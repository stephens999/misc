/* Stream a phased VCF on stdin and emit an LD-thinned, MAF-filtered set of
 * biallelic autosomal SNPs as derived-allele haplotype calls.
 *
 *   ./vcf_thin <min_dist_bp> <min_maf> <out_prefix>
 *
 * Writes <prefix>.samples (one sample id per line), <prefix>.snps
 * (chrom:pos ref alt aa daf) and <prefix>.hap, in which each line is one SNP
 * and holds 2*nsample characters '0'/'1', the derived-allele count of each
 * haplotype (sample 1 hap A, sample 1 hap B, sample 2 hap A, ...).
 *
 * Thinning is by physical distance: a variant is considered only if it is at
 * least min_dist from the last SNP kept on that chromosome, and the first
 * such variant that passes every filter is the one kept.  Filters are
 * FILTER == PASS, single-base REF and ALT (so biallelic SNPs only), an
 * AA_ensembl ancestral allele matching REF or ALT, no missing genotype in any
 * sample, and derived-allele frequency in [min_maf, 1 - min_maf].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define BUFSZ (1 << 22)

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: vcf_thin min_dist min_maf prefix\n"); return 1; }
    long  min_dist = atol(argv[1]);
    double min_maf = atof(argv[2]);
    const char *pre = argv[3];

    char path[4096];
    snprintf(path, sizeof path, "%s.samples", pre); FILE *fs = fopen(path, "w");
    snprintf(path, sizeof path, "%s.snps",    pre); FILE *fv = fopen(path, "w");
    snprintf(path, sizeof path, "%s.hap",     pre); FILE *fh = fopen(path, "w");
    if (!fs || !fv || !fh) { perror("fopen"); return 1; }

    char *line = malloc(BUFSZ);
    size_t cap = BUFSZ;
    ssize_t len;
    int nsamp = 0;
    char *hap = NULL;
    char last_chrom[64] = "";
    long last_pos = -1000000000L;
    long long nread = 0, nkept = 0;

    while ((len = getline(&line, &cap, stdin)) > 0) {
        if (line[0] == '#') {
            if (strncmp(line, "#CHROM", 6) == 0) {
                /* sample ids are fields 10.. */
                char *p = line; int f = 0;
                while (*p && *p != '\n') {
                    char *q = p; while (*q && *q != '\t' && *q != '\n') q++;
                    if (f >= 9) { fwrite(p, 1, q - p, fs); fputc('\n', fs); nsamp++; }
                    f++; p = (*q == '\t') ? q + 1 : q;
                }
                hap = malloc(2 * nsamp + 2);
                fprintf(stderr, "nsamples = %d\n", nsamp);
            }
            continue;
        }
        if (!nsamp) continue;
        nread++;

        /* CHROM */
        char *p = line;
        char *tab = strchr(p, '\t'); if (!tab) continue;
        size_t clen = tab - p; if (clen >= sizeof last_chrom) continue;
        char chrom[64]; memcpy(chrom, p, clen); chrom[clen] = 0;

        /* POS */
        p = tab + 1; long pos = strtol(p, &tab, 10); if (*tab != '\t') continue;

        if (strcmp(chrom, last_chrom) != 0) { strcpy(last_chrom, chrom); last_pos = -1000000000L; }
        if (pos - last_pos < min_dist) continue;          /* thinned out */

        /* ID */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        /* REF */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        if (tab - p != 1) continue;
        char ref = toupper(*p);
        /* ALT */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        if (tab - p != 1) continue;                        /* drops multiallelic */
        char alt = toupper(*p);
        if (strchr("ACGT", ref) == NULL || strchr("ACGT", alt) == NULL) continue;
        /* QUAL */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        /* FILTER */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        if (tab - p != 4 || strncmp(p, "PASS", 4) != 0) continue;
        /* INFO: ancestral allele */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        char *aap = strstr(p, "AA_ensembl=");
        if (!aap || aap > tab) continue;
        char aa = toupper(aap[11]);
        char aanext = aap[12];
        if (aanext != ';' && aanext != '\t' && aanext != '\n') continue;  /* multi-char */
        int derived_is_alt;
        if      (aa == ref) derived_is_alt = 1;
        else if (aa == alt) derived_is_alt = 0;
        else continue;                                      /* AA matches neither */
        /* FORMAT (must be plain GT) */
        p = tab + 1; tab = strchr(p, '\t'); if (!tab) continue;
        if (tab - p != 2 || strncmp(p, "GT", 2) != 0) continue;

        /* genotypes */
        p = tab + 1;
        int i, bad = 0; long nderiv = 0;
        for (i = 0; i < nsamp; i++) {
            char a = p[0], b = p[2];
            if (a == '.' || b == '.' || p[1] != '|') { bad = 1; break; }
            int da = (a == '1'), db = (b == '1');
            if (!derived_is_alt) { da = !da; db = !db; }
            hap[2 * i]     = '0' + da;
            hap[2 * i + 1] = '0' + db;
            nderiv += da + db;
            p += 3;
            if (i + 1 < nsamp) { if (*p != '\t') { bad = 1; break; } p++; }
        }
        if (bad) continue;

        double daf = (double) nderiv / (2.0 * nsamp);
        if (daf < min_maf || daf > 1.0 - min_maf) continue;

        last_pos = pos;
        nkept++;
        fprintf(fv, "%s:%ld\t%c\t%c\t%c\t%.6f\n", chrom, pos, ref, alt, aa, daf);
        fwrite(hap, 1, 2 * nsamp, fh); fputc('\n', fh);
        if (nkept % 10000 == 0) fprintf(stderr, "kept %lld of %lld at %s:%ld\n", nkept, nread, chrom, pos);
    }
    fprintf(stderr, "done: kept %lld of %lld records\n", nkept, nread);
    fclose(fs); fclose(fv); fclose(fh);
    return 0;
}
