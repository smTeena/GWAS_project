# Workflow of GWAS project

## 1. QUALITY CONTROL AND RELATEDNESS


### Filtering bfiles to only contain IDs with phenotypic data

```R
awk '{print $1, $1}' height.txt > keep_height_ids.txt
plink --bfile gwas_data --keep keep_height_ids.txt --make-bed --out gwas_data_height
#1376653 variants and 1071 people pass filters and QC.
```

### Filtering the height_metadata.txt and filtering by chip


```R
awk '{print $1}' height.txt > keep_height_ids_onecolumn.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids' keep_height_ids_onecolumn.txt metadata.txt > metadata_height.txt
wc -l metadata_height.txt
#1071 metadata_height.txt
```

```R
metadata <- read.table("metadata_height.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
```

 
  
    412   3671 37      OmniExpress v1 Female
    413   3934 37   HTS iSelect HD v4   Male
    414   3370 37                           
    415   3571 37      OmniExpress v1 Female
    416    347 37 OmniExpress plus v3 Female
    417    345 37 OmniExpress plus v3 Female
    418   3495 37   HTS iSelect HD v4 Female
 
    1060    Illumina GSAs
    1061    I1070      OniExpress
    1071      OmniExpress



```R
nrow(metadata)
```


1071



```R
for (chip_type in unique(metadata$chip)) {
  
  # Subset metadata by chip type
  group_data <- subset(metadata, chip == chip_type)
  
  # Clean chip_type string for filename
  filename <- paste0(gsub(" ", "_", chip_type), ".keep")
  
  # Write keep file (FID and IID as user ID)
  write.table(
    data.frame(FID = group_data$user, IID = group_data$user),
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
```

### Creating bfiles for each chip subset


```R
plink --bfile gwas_data_height \
      --keep HTS_iSelect_HD.keep \
      --make-bed \
      --out HTS_iSelect_HD
plink --bfile gwas_data_height \
      --keep Illumina_GSAs.keep \
      --make-bed \
      --out Illumina_GSAs
plink --bfile gwas_data_height \
      --keep Unknown.keep \
      --make-bed \
      --out Unknown
plink --bfile gwas_data_height \
      --keep OmniExpress_plus.keep \
      --make-bed \
      --out OmniExpress_plus
plink --bfile gwas_data_height \
      --keep OmniExpress.keep \
      --make-bed \
      --out OmniExpress
```

### Sample QC

#### 1. Sex check - removing all mismatches and unknowns


```R
plink --bfile HTS_iSelect_HD --check-sex --out HTS_iSelect_HD_sex
plink --bfile Illumina_GSAs --check-sex --out Illumina_GSAs_sex
plink --bfile Unknown --check-sex --out Unknown_sex
plink --bfile OmniExpress_plus --check-sex --out OmniExpress_plus_sex
plink --bfile OmniExpress --check-sex --out OmniExpress_sex

grep PROBLEM HTS_iSelect_HD_sex.sexcheck > HTS_iSelect_HD_wrong_sex.txt
plink --bfile HTS_iSelect_HD --remove HTS_iSelect_HD_wrong_sex.txt --make-bed --out HTS_iSelect_HD_postsex

grep PROBLEM Illumina_GSAs_sex.sexcheck > Illumina_GSAs_wrong_sex.txt
plink --bfile Illumina_GSAs --remove Illumina_GSAs_wrong_sex.txt --make-bed --out Illumina_GSAs_postsex

grep PROBLEM Unknown_sex.sexcheck > Unknown_wrong_sex.txt
plink --bfile Unknown --remove Unknown_wrong_sex.txt --make-bed --out Unknown_postsex

grep PROBLEM OmniExpress_plus_sex.sexcheck > OmniExpress_plus_wrong_sex.txt
plink --bfile OmniExpress_plus --remove OmniExpress_plus_wrong_sex.txt --make-bed --out OmniExpress_plus_postsex

grep PROBLEM OmniExpress_sex.sexcheck > OmniExpress_wrong_sex.txt
plink --bfile OmniExpress --remove OmniExpress_wrong_sex.txt --make-bed --out OmniExpress_postsex
```

#### 2. Identification of individuals with elevated missing data rates or outlying heterozygosity rate


```R
plink --bfile OmniExpress_postsex --missing --out OmniExpress_postsex
plink --bfile OmniExpress_postsex --het --out OmniExpress_postsex 

plink --bfile OmniExpress_plus_postsex --missing --out OmniExpress_plus_postsex
plink --bfile OmniExpress_plus_postsex --het --out OmniExpress_plus_postsex

plink --bfile Unknown_postsex --missing --out Unknown_postsex
plink --bfile Unknown_postsex --het --out Unknown_postsex

plink --bfile Illumina_GSAs_postsex --missing --out Illumina_GSAs_postsex
plink --bfile Illumina_GSAs_postsex --het --out Illumina_GSAs_postsex

plink --bfile HTS_iSelect_HD_postsex --missing --out HTS_iSelect_HD_postsex
plink --bfile HTS_iSelect_HD_postsex --het --out HTS_iSelect_HD_postsex
```

.imiss - Shows how much data is missing per person (i.e., per sample)
.lmiss - Shows how much data is missing per SNP.
.het - assess heterozygosity per individual
| Column | Description                                             |
| ------ | ------------------------------------------------------- |
| FID    | Family ID                                               |
| IID    | Individual ID                                           |
| O(HOM) | Observed number of homozygous genotypes                 |
| E(HOM) | Expected number of homozygous genotypes                 |
| N(NM)  | Number of non-missing genotypes                         |
| F      | Inbreeding coefficient estimate (observed vs. expected) |

F = 0 → Normal heterozygosity
F > 0 → More homozygous than expected (possible inbreeding)
F < 0 → More heterozygous than expected (possible contamination)


```R
library(dplyr)

# Read the missingness and heterozygosity files
HTS_iSelect_HD_postsex_miss <- read.table("HTS_iSelect_HD_postsex.imiss", header = TRUE)
HTS_iSelect_HD_postsex_het <- read.table("HTS_iSelect_HD_postsex.het", header = TRUE)

# Merge both datasets by FID and IID
HTS_iSelect_HD_postsex <- inner_join(HTS_iSelect_HD_postsex_miss, HTS_iSelect_HD_postsex_het)

# Calculate observed heterozygosity rate
HTS_iSelect_HD_postsex$Het <- (HTS_iSelect_HD_postsex$N.NM. - HTS_iSelect_HD_postsex$O.HOM.) / HTS_iSelect_HD_postsex$N.NM.

plot(HTS_iSelect_HD_postsex$Het, HTS_iSelect_HD_postsex$F_MISS, xlab="Heterozygosity Rate", ylab="Missing Rate", main="Heterozygosity vs Missingness (HTS_iSelect_HD)")

# Add threshold line for missingness (e.g., 0.02)

# Define heterozygosity thresholds (e.g., +/- 3 SD from mean)
het_mean <- mean(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
het_sd <- sd(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
abline(v = het_mean + 3 * het_sd, col = "blue", lty = 2)  # vertical upper threshold
abline(v = het_mean - 3 * het_sd, col = "blue", lty = 2)  # vertical lower threshold
abline(h = 0.02, col = "red", lty = 2)  # horizontal line
```

    [1m[22mJoining with `by = join_by(FID, IID)`



    
![png](output_20_1.png)
    


Extremely high missingness, possibly due to presence of SNPs that were not accounted in the chip specific array.


```R
HTS_lmiss <- read.table("HTS_iSelect_HD_postsex.lmiss", header = TRUE)
HTS_lmiss$chip <- "HTS_iSelect_HD"

Illumina_GSAs_lmiss <- read.table("Illumina_GSAs_postsex.lmiss", header = TRUE)
Illumina_GSAs_lmiss$chip <- "Illumina_GSAs"

nochip_lmiss <- read.table("Unknown_postsex.lmiss", header = TRUE)
nochip_lmiss$chip <- "Unknown"

OmniExpress_plus_lmiss <- read.table("OmniExpress_plus_postsex.lmiss", header = TRUE)
OmniExpress_plus_lmiss$chip <- "OmniExpress_plus"

OmniExpress_lmiss <- read.table("OmniExpress_postsex.lmiss", header = TRUE)
OmniExpress_lmiss$chip <- "OmniExpress"

all_lmiss <- rbind(HTS_lmiss, Illumina_GSAs_lmiss, nochip_lmiss, OmniExpress_plus_lmiss, OmniExpress_lmiss)
```


    
![png](output_22_0.png)
    



```R
library(ggplot2)

ggplot(all_lmiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", alpha = 0.7, boundary = 0) +
  geom_vline(xintercept = 1.0, color = "red", linetype = "dashed", size = 1) +  # 100% missing
  facet_wrap(~ chip, scales = "free_y") +
  theme_minimal() +
  labs(title = "SNP Missingness Distribution by Chip Type",
       subtitle = "Red dashed line marks SNPs with 100% missingness",
       x = "Fraction Missing (F_MISS)",
       y = "Number of SNPs") +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10))

```

    Warning message:
    “[1m[22mUsing `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    [36mℹ[39m Please use `linewidth` instead.”



    
![png](output_23_1.png)
    



```R
library(dplyr)

all_lmiss %>%
  filter(F_MISS == 1) %>%
  group_by(chip) %>%
  summarise("100%_missing_SNPs" = n())
```


<table class="dataframe">
<caption>A tibble: 5 × 2</caption>
<thead>
	<tr><th scope=col>chip</th><th scope=col>100%_missing_SNPs</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>HTS_iSelect_HD  </td><td>380169</td></tr>
	<tr><td>Illumina_GSAs   </td><td>853856</td></tr>
	<tr><td>OmniExpress     </td><td>   597</td></tr>
	<tr><td>OmniExpress_plus</td><td>342534</td></tr>
	<tr><td>Unknown         </td><td>     1</td></tr>
</tbody>
</table>




```R
library(dplyr)

# Loop through each chip and write list of SNPs with 100% missingness
unique_chips <- unique(all_lmiss$chip)

for (chip_name in unique_chips) {
  snps_to_remove <- all_lmiss %>%
    filter(chip == chip_name, F_MISS == 1) %>%
    pull(SNP)
  
  # Save to file
  write.table(snps_to_remove,
              file = paste0("snps_to_remove_", gsub(" ", "_", chip_name), ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

```

Removing the SNPs and creating new bfiles with chip relavant SNPs


```R
plink --bfile HTS_iSelect_HD_postsex --exclude snps_to_remove_HTS_iSelect_HD.txt --make-bed --out HTS_iSelect_HD_postsex_specificSNPs
plink --bfile Illumina_GSAs_postsex --exclude snps_to_remove_Illumina_GSAs.txt --make-bed --out Illumina_GSAs_postsex_specificSNPs
plink --bfile OmniExpress_postsex --exclude snps_to_remove_OmniExpress.txt --make-bed --out OmniExpress_postsex_specificSNPs
plink --bfile OmniExpress_plus_postsex --exclude snps_to_remove_OmniExpress_plus.txt --make-bed --out OmniExpress_plus_postsex_specificSNPs
plink --bfile Unknown_postsex --exclude snps_to_remove_Unknown.txt --make-bed --out Unknown_postsex_specificSNPs
```

Redoing the sample missingnes and heterozygosity


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs --missing --out HTS_iSelect_HD_postsex_specificSNPs
plink --bfile Illumina_GSAs_postsex_specificSNPs --missing --out Illumina_GSAs_postsex_specificSNPs
plink --bfile OmniExpress_postsex_specificSNPs --missing --out OmniExpress_postsex_specificSNPs
plink --bfile OmniExpress_plus_postsex_specificSNPs --missing --out OmniExpress_plus_postsex_specificSNPs
plink --bfile Unknown_postsex_specificSNPs --missing --out Unknown_postsex_specificSNPs

plink --bfile HTS_iSelect_HD_postsex_specificSNPs --het --out HTS_iSelect_HD_postsex_specificSNPs 
plink --bfile Illumina_GSAs_postsex_specificSNPs --het --out Illumina_GSAs_postsex_specificSNPs 
plink --bfile OmniExpress_postsex_specificSNPs --het --out OmniExpress_postsex_specificSNPs 
plink --bfile OmniExpress_plus_postsex_specificSNPs --het --out OmniExpress_plus_postsex_specificSNPs
plink --bfile Unknown_postsex_specificSNPs --het --out Unknown_postsex_specificSNPs

```


    Error in parse(text = x, srcfile = src): <text>:1:15: unexpected symbol
    1: plink --bfile HTS_iSelect_HD_postsex_specificSNPs
                      ^
    Traceback:




```R
library(dplyr)

# Read the missingness and heterozygosity files
HTS_iSelect_HD_postsex_miss <- read.table("HTS_iSelect_HD_postsex_specificSNPs.imiss", header = TRUE)
HTS_iSelect_HD_postsex_het <- read.table("HTS_iSelect_HD_postsex_specificSNPs.het", header = TRUE)

# Merge both datasets by FID and IID
HTS_iSelect_HD_postsex <- inner_join(HTS_iSelect_HD_postsex_miss, HTS_iSelect_HD_postsex_het)

# Calculate observed heterozygosity rate
HTS_iSelect_HD_postsex$Het <- (HTS_iSelect_HD_postsex$N.NM. - HTS_iSelect_HD_postsex$O.HOM.) / HTS_iSelect_HD_postsex$N.NM.

plot(HTS_iSelect_HD_postsex$Het, HTS_iSelect_HD_postsex$F_MISS, xlab="Heterozygosity Rate", ylab="Missing Rate", main="Heterozygosity vs Missingness (HTS_iSelect_HD)")

# Add threshold line for missingness (e.g., 0.02)

# Define heterozygosity thresholds (e.g., +/- 3 SD from mean)
het_mean <- mean(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
het_sd <- sd(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
abline(v = het_mean + 3 * het_sd, col = "blue", lty = 2)  # vertical upper threshold
abline(v = het_mean - 3 * het_sd, col = "blue", lty = 2)  # vertical lower threshold
abline(h = 0.02, col = "red", lty = 2)  # horizontal line
```

    [1m[22mJoining with `by = join_by(FID, IID)`



    
![png](output_30_1.png)
    



```R
library(dplyr)
library(ggplot2)

# Function to read and process each chip dataset
process_chip_data <- function(miss_file, het_file, chip_name) {
  miss <- read.table(miss_file, header = TRUE)
  het <- read.table(het_file, header = TRUE)
  
  df <- inner_join(miss, het, by = c("FID", "IID"))
  df <- df %>% 
    mutate(Het = (N.NM. - O.HOM.) / N.NM.,
           chip = chip_name)
  return(df)
}

# Process all chips (update file names accordingly)
chip1 <- process_chip_data("HTS_iSelect_HD_postsex_specificSNPs.imiss", "HTS_iSelect_HD_postsex_specificSNPs.het", "HTS iSelect HD")
chip2 <- process_chip_data("Illumina_GSAs_postsex_specificSNPs.imiss", "Illumina_GSAs_postsex_specificSNPs.het", "Illumina GSAs")
chip3 <- process_chip_data("OmniExpress_postsex_specificSNPs.imiss", "OmniExpress_postsex_specificSNPs.het", "OmniExpress")
chip4 <- process_chip_data("OmniExpress_plus_postsex_specificSNPs.imiss", "OmniExpress_plus_postsex_specificSNPs.het", "OmniExpress plus")
chip5 <- process_chip_data("Unknown_postsex_specificSNPs.imiss", "Unknown_postsex_specificSNPs.het", "Unknown")

# Combine all dataframes
all_data <- bind_rows(chip1, chip2, chip3, chip4, chip5)

# Calculate mean and SD per chip for heterozygosity thresholds
het_stats <- all_data %>%
  group_by(chip) %>%
  summarise(
    het_mean = mean(Het, na.rm = TRUE),
    het_sd = sd(Het, na.rm = TRUE)
  )

# Join thresholds back to the main dataframe
all_data <- all_data %>%
  left_join(het_stats, by = "chip")

# Plot with ggplot2 and facet_wrap
ggplot(all_data, aes(x = Het, y = F_MISS)) +
  geom_point(alpha = 0.6) +
  geom_vline(aes(xintercept = het_mean + 3 * het_sd), color = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = het_mean - 3 * het_sd), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed") +
  facet_wrap(~ chip, scales = "free") +
  labs(title = "Heterozygosity vs Missingness by Chip",
       subtitle = "After removing SNPs with 100% missiness",
       x = "Heterozygosity Rate",
       y = "Missing Rate") +
  theme_minimal()

```

    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](output_31_1.png)
    


Missingness is still bad


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs --geno 0.05 --make-bed --out HTS_iSelect_HD_postsex_specificSNPs_geno0.05
plink --bfile Illumina_GSAs_postsex_specificSNPs --geno 0.05 --make-bed --out Illumina_GSAs_postsex_specificSNPs_geno0.05
plink --bfile OmniExpress_plus_postsex_specificSNPs --geno 0.05 --make-bed --out OmniExpress_plus_postsex_specificSNPs_geno0.05
plink --bfile OmniExpress_postsex_specificSNPs --geno 0.05 --make-bed --out OmniExpress_postsex_specificSNPs_geno0.05
plink --bfile Unknown_postsex_specificSNPs --geno 0.05 --make-bed --out Unknown_postsex_specificSNPs_geno0.05
```


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno0.05 --missing --out HTS_iSelect_HD_postsex_specificSNPs_geno0.05
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno0.05 --het --out HTS_iSelect_HD_postsex_specificSNPs_geno0.05 

plink --bfile Illumina_GSAs_postsex_specificSNPs_geno0.05 --missing --out Illumina_GSAs_postsex_specificSNPs_geno0.05
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno0.05 --het --out Illumina_GSAs_postsex_specificSNPs_geno0.05 

plink --bfile OmniExpress_plus_postsex_specificSNPs_geno0.05 --missing --out OmniExpress_plus_postsex_specificSNPs_geno0.05
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno0.05 --het --out OmniExpress_plus_postsex_specificSNPs_geno0.05 

plink --bfile OmniExpress_postsex_specificSNPs_geno0.05 --missing --out OmniExpress_postsex_specificSNPs_geno0.05
plink --bfile OmniExpress_postsex_specificSNPs_geno0.05 --het --out OmniExpress_postsex_specificSNPs_geno0.05

plink --bfile Unknown_postsex_specificSNPs_geno0.05 --missing --out Unknown_postsex_specificSNPs_geno0.05
plink --bfile Unknown_postsex_specificSNPs_geno0.05 --het --out Unknown_postsex_specificSNPs_geno0.05
```


```R
library(dplyr)

# Read the missingness and heterozygosity files
HTS_iSelect_HD_postsex_miss <- read.table("HTS_iSelect_HD_postsex_specificSNPs_geno0.05.imiss", header = TRUE)
HTS_iSelect_HD_postsex_het <- read.table("HTS_iSelect_HD_postsex_specificSNPs_geno0.05.het", header = TRUE)

# Merge both datasets by FID and IID
HTS_iSelect_HD_postsex <- inner_join(HTS_iSelect_HD_postsex_miss, HTS_iSelect_HD_postsex_het)

# Calculate observed heterozygosity rate
HTS_iSelect_HD_postsex$Het <- (HTS_iSelect_HD_postsex$N.NM. - HTS_iSelect_HD_postsex$O.HOM.) / HTS_iSelect_HD_postsex$N.NM.

plot(HTS_iSelect_HD_postsex$Het, HTS_iSelect_HD_postsex$F_MISS, xlab="Heterozygosity Rate", ylab="Missing Rate", main="Heterozygosity vs Missingness (HTS_iSelect_HD)")

# Add threshold line for missingness (e.g., 0.02)

# Define heterozygosity thresholds (e.g., +/- 3 SD from mean)
het_mean <- mean(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
het_sd <- sd(HTS_iSelect_HD_postsex$Het, na.rm = TRUE)
abline(v = het_mean + 3 * het_sd, col = "blue", lty = 2)  # vertical upper threshold
abline(v = het_mean - 3 * het_sd, col = "blue", lty = 2)  # vertical lower threshold
abline(h = 0.02, col = "red", lty = 2)  # horizontal line
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    [1m[22mJoining with `by = join_by(FID, IID)`



    
![png](output_35_1.png)
    



```R
library(dplyr)
library(ggplot2)

# Function to load and process each chip-specific dataset
process_chip_data <- function(miss_file, het_file, chip_name) {
  miss <- read.table(miss_file, header = TRUE)
  het <- read.table(het_file, header = TRUE)
  
  df <- inner_join(miss, het, by = c("FID", "IID")) %>%
    mutate(
      Het = (N.NM. - O.HOM.) / N.NM.,
      chip = chip_name
    )
  
  return(df)
}

# Load all datasets (edit file paths and chip names accordingly)
chip1 <- process_chip_data("HTS_iSelect_HD_postsex_specificSNPs_geno0.05.imiss", "HTS_iSelect_HD_postsex_specificSNPs_geno0.05.het", "HTS iSelect HD")
chip2 <- process_chip_data("Illumina_GSAs_postsex_specificSNPs_geno0.05.imiss", "Illumina_GSAs_postsex_specificSNPs_geno0.05.het", "Illumina GSAs")
chip3 <- process_chip_data("OmniExpress_postsex_specificSNPs_geno0.05.imiss", "OmniExpress_postsex_specificSNPs_geno0.05.het", "OmniExpress")
chip4 <- process_chip_data("OmniExpress_plus_postsex_specificSNPs_geno0.05.imiss", "OmniExpress_plus_postsex_specificSNPs_geno0.05.het", "OmniExpress plus")
chip5 <- process_chip_data("Unknown_postsex_specificSNPs_geno0.05.imiss", "Unknown_postsex_specificSNPs_geno0.05.het", "Unknown")

# Combine datasets
all_data <- bind_rows(chip1, chip2, chip3, chip4, chip5)

# Compute mean and SD of heterozygosity per chip
het_stats <- all_data %>%
  group_by(chip) %>%
  summarise(
    het_mean = mean(Het, na.rm = TRUE),
    het_sd = sd(Het, na.rm = TRUE)
  )

# Merge stats back to main dataframe
all_data <- left_join(all_data, het_stats, by = "chip")

# Plot with ggplot2
ggplot(all_data, aes(x = Het, y = F_MISS)) +
  geom_point(alpha = 0.6) +
  geom_vline(aes(xintercept = het_mean + 3 * het_sd), color = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = het_mean - 3 * het_sd), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 0.02, color = "red", linetype = "dashed") +
  facet_wrap(~ chip, scales = "free") +
  theme_minimal() +
  labs(
    title = "Heterozygosity vs Missingness (post --geno 0.05)",
    x = "Heterozygosity Rate",
    y = "Missing Genotype Rate"
  )

```

    Warning message:
    “[1m[22mRemoved 6 rows containing missing values (`geom_point()`).”



    
![png](output_36_1.png)
    


Remove heterozygosity +- 3 s.d. from mean


```R
# Calculate upper and lower Het thresholds
upper_het_thresh <- het_mean + 3 * het_sd
lower_het_thresh <- het_mean - 3 * het_sd

# Identify samples to exclude based on Het
het_outliers <- HTS_iSelect_HD_postsex %>%
  filter(Het > upper_het_thresh | Het < lower_het_thresh) %>%
  select(FID, IID)

# Write out a file listing these samples for removal
write.table(het_outliers, "HTS_iSelect_HD_het_outliers.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
```


```R
library(dplyr)

# Define a function for heterozygosity filtering
process_het_outliers <- function(miss_file, het_file, output_name) {
  # Load files
  miss <- read.table(miss_file, header = TRUE)
  het <- read.table(het_file, header = TRUE)

  # Merge and compute heterozygosity
  merged <- inner_join(miss, het, by = c("FID", "IID")) %>%
    mutate(Het = (N.NM. - O.HOM.) / N.NM.)

  # Calculate thresholds
  het_mean <- mean(merged$Het, na.rm = TRUE)
  het_sd <- sd(merged$Het, na.rm = TRUE)
  upper_het_thresh <- het_mean + 3 * het_sd
  lower_het_thresh <- het_mean - 3 * het_sd

  # Filter outliers
  het_outliers <- merged %>%
    filter(Het > upper_het_thresh | Het < lower_het_thresh) %>%
    select(FID, IID)

  # Save to file
  write.table(het_outliers, paste0(output_name, "_het_outliers.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Return count of outliers (optional)
  return(nrow(het_outliers))
}

# Illumina GSAs
process_het_outliers("Illumina_GSAs_postsex_specificSNPs_geno0.05.imiss",
                     "Illumina_GSAs_postsex_specificSNPs_geno0.05.het",
                     "Illumina_GSAs")

# OmniExpress
process_het_outliers("OmniExpress_postsex_specificSNPs_geno0.05.imiss",
                     "OmniExpress_postsex_specificSNPs_geno0.05.het",
                     "OmniExpress")

# OmniExpress Plus
process_het_outliers("OmniExpress_plus_postsex_specificSNPs_geno0.05.imiss",
                     "OmniExpress_plus_postsex_specificSNPs_geno0.05.het",
                     "OmniExpress_plus")
#Unknown
process_het_outliers("Unknown_postsex_specificSNPs_geno0.05.imiss",
                     "Unknown_postsex_specificSNPs_geno0.05.het",
                     "Unknown")
```


1



4



1



1



```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno0.05 --remove HTS_iSelect_HD_het_outliers.txt --make-bed --out HTS_iSelect_HD_postsex_specificSNPs_geno_het
#262 people remaining
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno0.05 --remove Illumina_GSAs_het_outliers.txt --make-bed --out Illumina_GSAs_postsex_specificSNPs_geno_het
#124 people remaining.
plink --bfile OmniExpress_postsex_specificSNPs_geno0.05 --remove OmniExpress_het_outliers.txt --make-bed --out OmniExpress_postsex_specificSNPs_geno_het
#125 people remaining.
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno0.05 --remove OmniExpress_plus_het_outliers.txt --make-bed --out OmniExpress_plus_postsex_specificSNPs_geno_het
#257 people remaining.
plink --bfile Unknown_postsex_specificSNPs_geno0.05 --remove Unknown_het_outliers.txt --make-bed --out Unknown_postsex_specificSNPs_geno_het
#202 people remaining.
```


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het --mind 0.02 --make-bed --out HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind
#260 people remaining
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het --mind 0.02 --make-bed --out Illumina_GSAs_postsex_specificSNPs_geno_het_mind
#124 people remaining.
plink --bfile OmniExpress_postsex_specificSNPs_geno_het --mind 0.02 --make-bed --out OmniExpress_postsex_specificSNPs_geno_het_mind
#121 people remaining.
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het --mind 0.02 --make-bed --out OmniExpress_plus_postsex_specificSNPs_geno_het_mind
#256 people remaining.
plink --bfile Unknown_postsex_specificSNPs_geno_het --mind 0.02 --make-bed --out Unknown_postsex_specificSNPs_geno_het_mind
#176 people remaining.
```

### Identification of duplicated or related individuals

LD prune the SNPs to get a list of independent SNPs to keep


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind --indep-pairwise 500kb 5 0.2 --out HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind
#429870 of 571075 variants removed
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het_mind --indep-pairwise 500kb 5 0.2 --out Illumina_GSAs_postsex_specificSNPs_geno_het_mind
#323678 of 514384 variants removed.
plink --bfile OmniExpress_postsex_specificSNPs_geno_het_mind --indep-pairwise 500kb 5 0.2 --out OmniExpress_postsex_specificSNPs_geno_het_mind
#521916 of 638800 variants removed.
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het_mind --indep-pairwise 500kb 5 0.2 --out OmniExpress_plus_postsex_specificSNPs_geno_het_mind
#773657 of 913802 variants removed.
plink --bfile Unknown_postsex_specificSNPs_geno_het_mind --indep-pairwise 500kb 5 0.2 --out Unknown_postsex_specificSNPs_geno_het_mind
#23702 of 66864 variants removed.

```


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind --extract HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind.prune.in --genome --min 0.185 --out HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_ldprune
# 0 related
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het_mind --extract Illumina_GSAs_postsex_specificSNPs_geno_het_mind.prune.in --genome --min 0.185 --out Illumina_GSAs_postsex_specificSNPs_geno_het_mind_ldprune
# 1 pair related
plink --bfile OmniExpress_postsex_specificSNPs_geno_het_mind --extract OmniExpress_postsex_specificSNPs_geno_het_mind.prune.in --genome --min 0.185 --out OmniExpress_postsex_specificSNPs_geno_het_mind_ldprune
# 2 pairs related
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het_mind --extract OmniExpress_plus_postsex_specificSNPs_geno_het_mind.prune.in --genome --min 0.185 --out OmniExpress_plus_postsex_specificSNPs_geno_het_mind_ldprune
# 1 pair related
plink --bfile Unknown_postsex_specificSNPs_geno_het_mind --extract Unknown_postsex_specificSNPs_geno_het_mind.prune.in --genome --min 0.185 --out Unknown_postsex_specificSNPs_geno_het_mind_ldprune
# 2 pairs related
```


```R
# Manually changed names since function approached resulted in encoding error
# Set working directory or full paths
genome_file <- "OmniExpress_plus_postsex_specificSNPs_geno_het_mind_ldprune.genome"

# Read the .genome file
related <- read.table(genome_file, header = TRUE)

# Filter for related pairs (PI_HAT > 0.185)
related_pairs <- related[related$PI_HAT > 0.185, c("FID1", "IID1", "FID2", "IID2")]

# Choose one individual per pair to remove (e.g., always IID2)
to_remove <- unique(data.frame(FID = related_pairs$FID2, IID = related_pairs$IID2))

# Write to file for PLINK removal (UTF-8 safe)
write.table(to_remove, file = "OmniExpress_plus_related_pairs_to_remove.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding = "UTF-8")

```


```R

```

combine all the related info


```R
library(dplyr)

# List your .genome files (adjust names as needed)
genome_files <- c(
    "Unknown_postsex_specificSNPs_geno_het_mind_ldprune.genome",
    "OmniExpress_plus_postsex_specificSNPs_geno_het_mind_ldprune.genome",
    "OmniExpress_postsex_specificSNPs_geno_het_mind_ldprune.genome",
    "Illumina_GSAs_postsex_specificSNPs_geno_het_mind_ldprune.genome",
    "HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_ldprune.genome"
)

# Read and combine them
all_genome <- do.call(rbind, lapply(genome_files, function(file) {
  read.table(file, header = TRUE)
}))

# Optionally, write out the combined file
write.table(all_genome, "combined.genome", row.names = FALSE, quote = FALSE)

```

Removing one individual from the related pair and new bfiles


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind --remove HTS_iSelect_HD_related_pairs_to_remove.txt --make-bed --out HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_unrelated
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het_mind --remove OmniExpress_plus_related_pairs_to_remove.txt --make-bed --out OmniExpress_plus_postsex_specificSNPs_geno_het_mind_unrelated
plink --bfile OmniExpress_postsex_specificSNPs_geno_het_mind --remove OmniExpress_related_pairs_to_remove.txt --make-bed --out OmniExpress_postsex_specificSNPs_geno_het_mind_unrelated
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het_mind --remove Illumina_GSAs_related_pairs_to_remove.txt --make-bed --out Illumina_GSAs_postsex_specificSNPs_geno_het_mind_unrelated
plink --bfile Unknown_postsex_specificSNPs_geno_het_mind --remove Unknown_related_pairs_to_remove.txt --make-bed --out Unknown_postsex_specificSNPs_geno_het_mind_unrelated
```

SNP QC (hwe and maf)


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_unrelated --hwe 0.00001 --maf 0.01 --make-bed --out HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het_mind_unrelated --hwe 0.00001 --maf 0.01 --make-bed --out OmniExpress_plus_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf
plink --bfile OmniExpress_postsex_specificSNPs_geno_het_mind_unrelated --hwe 0.00001 --maf 0.01 --make-bed --out OmniExpress_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het_mind_unrelated --hwe 0.00001 --maf 0.01 --make-bed --out Illumina_GSAs_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf
plink --bfile Unknown_postsex_specificSNPs_geno_het_mind_unrelated --hwe 0.00001 --maf 0.01 --make-bed --out Unknown_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf
```

merge different chips


```R
# Load required library
library(dplyr)

# Read SNP (column 2) from each .bim file
HTS_iSelect_HD <- read.table("HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf.bim", stringsAsFactors = FALSE)[,2]
OmniExpress_plus <- read.table("OmniExpress_plus_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf.bim", stringsAsFactors = FALSE)[,2]
OmniExpress <- read.table("OmniExpress_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf.bim", stringsAsFactors = FALSE)[,2]
Illumina_GSAs <- read.table("Illumina_GSAs_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf.bim", stringsAsFactors = FALSE)[,2]
Unknown <- read.table("Unknown_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf.bim", stringsAsFactors = FALSE)[,2]

# Find intersection
common_snps <- Reduce(intersect, list(HTS_iSelect_HD, OmniExpress_plus, OmniExpress, Illumina_GSAs, Unknown))

# Save to file
write.table(common_snps, "common_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

```


```R
plink --bfile HTS_iSelect_HD_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf \
      --extract common_snps.txt \
      --make-bed \
      --out HTS_iSelect_HD_common
plink --bfile OmniExpress_plus_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf \
      --extract common_snps.txt \
      --make-bed \
      --out OmniExpress_plus_common
plink --bfile OmniExpress_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf \
      --extract common_snps.txt \
      --make-bed \
      --out OmniExpress_common
plink --bfile Illumina_GSAs_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf \
      --extract common_snps.txt \
      --make-bed \
      --out Illumina_GSAs_common
plink --bfile Unknown_postsex_specificSNPs_geno_het_mind_unrelated_hwemaf \
      --extract common_snps.txt \
      --make-bed \
      --out Unknown_common
```


```R
#created a merge_list.txt with all chips' bfiles listed
plink --bfile HTS_iSelect_HD_common \
      --merge-list merge_list.txt \
      --make-bed \
      --out merged_allchips_common

```

get PCs and plot for PCA


```R
plink --bfile merged_allchips_common --pca 20 --out merged_allchips_common_pca
```


```R
pcs <- read.table('merged_allchips_common_pca.eigenvec', head=F)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:(ncol(pcs)-2)))
pcs
```


<table class="dataframe">
<caption>A data.frame: 931 × 22</caption>
<thead>
	<tr><th scope=col>FID</th><th scope=col>IID</th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>⋯</th><th scope=col>PC11</th><th scope=col>PC12</th><th scope=col>PC13</th><th scope=col>PC14</th><th scope=col>PC15</th><th scope=col>PC16</th><th scope=col>PC17</th><th scope=col>PC18</th><th scope=col>PC19</th><th scope=col>PC20</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>  8</td><td>  8</td><td>-0.00695248</td><td>-0.000978118</td><td>-0.01889710</td><td>-0.040866000</td><td>-0.002266610</td><td> 0.03517340</td><td>-0.000920792</td><td>-0.051659200</td><td>⋯</td><td> 0.003243350</td><td> 0.01337570</td><td>-0.01707750</td><td> 0.007659100</td><td> 0.016004200</td><td> 0.027553000</td><td>-0.018682500</td><td> 0.01440690</td><td>-0.007379570</td><td> 0.02512950</td></tr>
	<tr><td> 11</td><td> 11</td><td> 0.00638909</td><td>-0.000885972</td><td>-0.06231530</td><td>-0.114144000</td><td>-0.023210300</td><td> 0.06640400</td><td>-0.011471100</td><td>-0.050189500</td><td>⋯</td><td>-0.045475600</td><td> 0.03473700</td><td> 0.00277825</td><td> 0.032124800</td><td>-0.044188900</td><td> 0.006255820</td><td>-0.014364400</td><td> 0.01456890</td><td> 0.017087500</td><td> 0.02230760</td></tr>
	<tr><td> 14</td><td> 14</td><td>-0.01231010</td><td>-0.004080500</td><td> 0.00722458</td><td> 0.008922570</td><td> 0.004325150</td><td>-0.02162320</td><td> 0.027711000</td><td> 0.042151200</td><td>⋯</td><td>-0.008133420</td><td> 0.02157040</td><td> 0.01742410</td><td> 0.023289400</td><td>-0.035920200</td><td> 0.006674540</td><td>-0.017092000</td><td>-0.02619020</td><td>-0.014704500</td><td>-0.01178390</td></tr>
	<tr><td> 16</td><td> 16</td><td>-0.01386140</td><td>-0.002910680</td><td>-0.00355122</td><td> 0.014365100</td><td>-0.011124600</td><td> 0.03915320</td><td>-0.002504390</td><td>-0.021549200</td><td>⋯</td><td>-0.005304090</td><td> 0.01849060</td><td>-0.02423420</td><td> 0.018379900</td><td> 0.000434086</td><td>-0.076599200</td><td>-0.017849000</td><td>-0.00166217</td><td> 0.004075230</td><td> 0.02687590</td></tr>
	<tr><td> 17</td><td> 17</td><td>-0.01186810</td><td>-0.001851780</td><td> 0.01857800</td><td> 0.027578300</td><td>-0.001946620</td><td>-0.02172390</td><td>-0.013766200</td><td> 0.039818800</td><td>⋯</td><td>-0.012517700</td><td>-0.00879131</td><td> 0.01622040</td><td>-0.002328350</td><td> 0.017400600</td><td>-0.032428600</td><td> 0.044414400</td><td> 0.00827015</td><td> 0.032261000</td><td>-0.01576650</td></tr>
	<tr><td> 26</td><td> 26</td><td>-0.00637205</td><td> 0.002338360</td><td> 0.01158400</td><td> 0.018724600</td><td> 0.018373300</td><td>-0.02935470</td><td>-0.005018150</td><td>-0.016958800</td><td>⋯</td><td>-0.002515420</td><td>-0.01636600</td><td> 0.01959510</td><td>-0.008870700</td><td>-0.018862000</td><td> 0.033481800</td><td> 0.006790150</td><td>-0.02563430</td><td> 0.021596300</td><td> 0.00792364</td></tr>
	<tr><td> 33</td><td> 33</td><td> 0.06932730</td><td> 0.115638000</td><td>-0.13813100</td><td> 0.093465900</td><td> 0.013792200</td><td>-0.00352651</td><td> 0.018360200</td><td>-0.016301000</td><td>⋯</td><td>-0.048628900</td><td> 0.01839450</td><td>-0.02686200</td><td>-0.028025900</td><td>-0.014875200</td><td>-0.146786000</td><td> 0.018922700</td><td>-0.03233660</td><td> 0.002884610</td><td> 0.00295314</td></tr>
	<tr><td> 35</td><td> 35</td><td>-0.01021090</td><td>-0.004794650</td><td>-0.00181076</td><td>-0.001693630</td><td>-0.015542000</td><td> 0.04344090</td><td> 0.022303000</td><td>-0.023907000</td><td>⋯</td><td> 0.024084400</td><td>-0.04464410</td><td>-0.01675030</td><td>-0.037671800</td><td> 0.002437710</td><td> 0.037906500</td><td> 0.016100800</td><td>-0.01066880</td><td> 0.000521250</td><td>-0.00637508</td></tr>
	<tr><td> 40</td><td> 40</td><td> 0.01429730</td><td> 0.003733020</td><td> 0.03214610</td><td>-0.025144700</td><td> 0.056398100</td><td>-0.00105906</td><td> 0.021974400</td><td> 0.007635840</td><td>⋯</td><td>-0.067036200</td><td>-0.01060720</td><td> 0.05504780</td><td>-0.079483800</td><td> 0.007869250</td><td> 0.013200800</td><td>-0.040183000</td><td>-0.03665850</td><td> 0.027211600</td><td> 0.02687350</td></tr>
	<tr><td> 60</td><td> 60</td><td>-0.00378555</td><td> 0.000465776</td><td>-0.01809920</td><td>-0.046011400</td><td>-0.002646640</td><td> 0.00812926</td><td>-0.009401230</td><td> 0.003673210</td><td>⋯</td><td>-0.052433500</td><td>-0.00607616</td><td> 0.02373270</td><td> 0.009950220</td><td>-0.036342900</td><td> 0.001114410</td><td> 0.033286000</td><td>-0.00133320</td><td> 0.010591200</td><td>-0.04244950</td></tr>
	<tr><td> 64</td><td> 64</td><td>-0.00343435</td><td>-0.005524260</td><td>-0.02668920</td><td>-0.062898200</td><td>-0.016058400</td><td> 0.08548320</td><td>-0.034246900</td><td> 0.013175100</td><td>⋯</td><td>-0.085484500</td><td>-0.04066220</td><td> 0.01374600</td><td> 0.006454250</td><td> 0.000728047</td><td>-0.024132800</td><td>-0.040959200</td><td> 0.02833070</td><td> 0.043158500</td><td> 0.02029640</td></tr>
	<tr><td> 71</td><td> 71</td><td>-0.01215810</td><td>-0.003156830</td><td> 0.02109350</td><td> 0.009314080</td><td> 0.000470235</td><td> 0.00361880</td><td> 0.020903100</td><td>-0.007933030</td><td>⋯</td><td> 0.002241200</td><td>-0.01639240</td><td>-0.00152056</td><td>-0.023761600</td><td> 0.009225780</td><td>-0.054128600</td><td> 0.018075600</td><td>-0.03406010</td><td> 0.005577290</td><td>-0.00206169</td></tr>
	<tr><td> 80</td><td> 80</td><td>-0.00813685</td><td>-0.004205830</td><td>-0.00643029</td><td>-0.012908400</td><td> 0.007847490</td><td>-0.01392740</td><td>-0.050897000</td><td> 0.032593400</td><td>⋯</td><td>-0.028844000</td><td>-0.01743010</td><td> 0.06745010</td><td> 0.010411000</td><td> 0.035194000</td><td>-0.000181053</td><td> 0.006786940</td><td> 0.00086782</td><td> 0.037219700</td><td>-0.02137390</td></tr>
	<tr><td> 86</td><td> 86</td><td> 0.06612750</td><td> 0.102072000</td><td>-0.19247100</td><td> 0.088194500</td><td> 0.047903300</td><td>-0.01776470</td><td> 0.005877560</td><td>-0.000926718</td><td>⋯</td><td> 0.036086700</td><td> 0.04015850</td><td> 0.05126590</td><td>-0.101954000</td><td> 0.005391420</td><td> 0.081873700</td><td> 0.021558200</td><td> 0.02484700</td><td>-0.005198440</td><td>-0.00791060</td></tr>
	<tr><td> 99</td><td> 99</td><td>-0.01012290</td><td>-0.008428760</td><td> 0.00749953</td><td> 0.004204310</td><td>-0.013123800</td><td>-0.00315479</td><td> 0.048796300</td><td>-0.007844710</td><td>⋯</td><td>-0.002723130</td><td> 0.02856540</td><td>-0.02500000</td><td>-0.025397400</td><td>-0.027949300</td><td>-0.013087900</td><td> 0.069385400</td><td> 0.02358360</td><td> 0.017773900</td><td>-0.00726219</td></tr>
	<tr><td>123</td><td>123</td><td>-0.01289580</td><td>-0.006148440</td><td> 0.00798913</td><td> 0.017812600</td><td>-0.001093100</td><td>-0.01874620</td><td> 0.003648930</td><td> 0.026978200</td><td>⋯</td><td>-0.010144800</td><td>-0.03321980</td><td> 0.02046720</td><td> 0.021591500</td><td> 0.040639100</td><td>-0.003499400</td><td>-0.011247300</td><td> 0.01219680</td><td>-0.028276400</td><td> 0.02903980</td></tr>
	<tr><td>137</td><td>137</td><td>-0.00650914</td><td> 0.001267750</td><td> 0.00233278</td><td>-0.010151300</td><td> 0.000740280</td><td>-0.02606450</td><td>-0.007467220</td><td>-0.079673800</td><td>⋯</td><td> 0.002548310</td><td> 0.02728730</td><td>-0.01491840</td><td>-0.022100300</td><td>-0.018257300</td><td>-0.056697100</td><td>-0.011568100</td><td>-0.03050540</td><td> 0.043137200</td><td>-0.00406372</td></tr>
	<tr><td>139</td><td>139</td><td>-0.01061170</td><td>-0.008603500</td><td> 0.01167920</td><td>-0.000667131</td><td> 0.005667180</td><td>-0.01325090</td><td>-0.023779000</td><td> 0.001530580</td><td>⋯</td><td>-0.022336600</td><td> 0.01620480</td><td>-0.03378190</td><td> 0.017014700</td><td>-0.003992730</td><td> 0.026869300</td><td>-0.036284500</td><td>-0.02712030</td><td> 0.045011200</td><td> 0.04572240</td></tr>
	<tr><td>141</td><td>141</td><td>-0.01226570</td><td>-0.007570240</td><td> 0.01160390</td><td> 0.000836282</td><td> 0.003740950</td><td>-0.02618460</td><td> 0.008483680</td><td> 0.030259700</td><td>⋯</td><td>-0.013619800</td><td> 0.00902871</td><td> 0.01277550</td><td>-0.000938854</td><td>-0.008100140</td><td>-0.024340800</td><td>-0.060282600</td><td> 0.00901595</td><td>-0.017402000</td><td>-0.00410366</td></tr>
	<tr><td>146</td><td>146</td><td> 0.06456810</td><td> 0.127522000</td><td> 0.06545520</td><td>-0.013323600</td><td>-0.035023300</td><td>-0.00225025</td><td>-0.064405200</td><td>-0.050655200</td><td>⋯</td><td>-0.004214050</td><td> 0.01739970</td><td> 0.00378024</td><td> 0.038588000</td><td>-0.000705364</td><td> 0.159936000</td><td>-0.062483500</td><td> 0.01403170</td><td>-0.076321300</td><td> 0.02374280</td></tr>
	<tr><td>154</td><td>154</td><td> 0.05383790</td><td> 0.100130000</td><td> 0.10587700</td><td>-0.061981000</td><td> 0.291509000</td><td> 0.05080930</td><td> 0.046581900</td><td>-0.014846900</td><td>⋯</td><td> 0.028406900</td><td> 0.02560340</td><td> 0.03588500</td><td>-0.004022170</td><td>-0.129999000</td><td>-0.012566300</td><td> 0.000779795</td><td>-0.01587260</td><td>-0.071997700</td><td>-0.02774990</td></tr>
	<tr><td>159</td><td>159</td><td>-0.01306500</td><td>-0.003161380</td><td> 0.01293670</td><td> 0.019311300</td><td>-0.006486250</td><td>-0.03421450</td><td> 0.023330200</td><td> 0.029786500</td><td>⋯</td><td>-0.000492562</td><td>-0.02818240</td><td> 0.00417504</td><td> 0.022899000</td><td>-0.014208300</td><td> 0.021784200</td><td> 0.023660500</td><td>-0.00404878</td><td>-0.053892900</td><td> 0.03137810</td></tr>
	<tr><td>160</td><td>160</td><td> 0.14668300</td><td>-0.100057000</td><td> 0.01269950</td><td> 0.000447172</td><td> 0.016715200</td><td>-0.00859451</td><td> 0.065884400</td><td>-0.067829100</td><td>⋯</td><td> 0.036803700</td><td> 0.09012120</td><td> 0.09345230</td><td> 0.092184200</td><td> 0.033766600</td><td> 0.125802000</td><td>-0.071428100</td><td> 0.10578200</td><td> 0.082100400</td><td> 0.07256560</td></tr>
	<tr><td>177</td><td>177</td><td>-0.00890694</td><td> 0.000920759</td><td>-0.00988666</td><td> 0.001172310</td><td>-0.007225970</td><td> 0.00736285</td><td>-0.003517460</td><td> 0.032475800</td><td>⋯</td><td> 0.081048300</td><td> 0.01341540</td><td>-0.00241616</td><td> 0.026697600</td><td>-0.049473700</td><td>-0.022813400</td><td> 0.038479200</td><td> 0.02521900</td><td>-0.005670120</td><td> 0.01527560</td></tr>
	<tr><td>180</td><td>180</td><td>-0.01157560</td><td>-0.002965810</td><td> 0.00462113</td><td> 0.011012900</td><td>-0.004538330</td><td>-0.03783210</td><td> 0.027494100</td><td>-0.022652600</td><td>⋯</td><td> 0.027065000</td><td>-0.05324650</td><td> 0.02668810</td><td> 0.031632000</td><td>-0.003453930</td><td> 0.025256700</td><td>-0.041711000</td><td>-0.04789800</td><td>-0.028999800</td><td> 0.01251820</td></tr>
	<tr><td>187</td><td>187</td><td>-0.00887721</td><td>-0.001800330</td><td> 0.01460560</td><td> 0.047367800</td><td>-0.013293100</td><td> 0.05935010</td><td>-0.039252300</td><td> 0.031670500</td><td>⋯</td><td> 0.017721200</td><td>-0.00818652</td><td> 0.00995422</td><td> 0.006483420</td><td>-0.011002200</td><td> 0.042867700</td><td> 0.005365060</td><td> 0.05113070</td><td>-0.006598700</td><td> 0.03147980</td></tr>
	<tr><td>199</td><td>199</td><td>-0.01301180</td><td>-0.006598700</td><td> 0.01271570</td><td> 0.018840400</td><td> 0.002226100</td><td>-0.01440890</td><td> 0.028451200</td><td> 0.018003500</td><td>⋯</td><td> 0.014580500</td><td> 0.03377430</td><td>-0.02987440</td><td>-0.010621100</td><td>-0.028464100</td><td>-0.017247000</td><td>-0.026144600</td><td> 0.01824580</td><td> 0.008674910</td><td> 0.01437120</td></tr>
	<tr><td>216</td><td>216</td><td>-0.01134360</td><td>-0.007507630</td><td> 0.00954351</td><td> 0.026520000</td><td> 0.002146290</td><td> 0.09247040</td><td> 0.030095400</td><td> 0.020834000</td><td>⋯</td><td> 0.017820000</td><td> 0.01926220</td><td>-0.01853720</td><td> 0.008214080</td><td>-0.020918600</td><td>-0.016168000</td><td>-0.095454200</td><td> 0.01688610</td><td> 0.000411898</td><td> 0.02518420</td></tr>
	<tr><td>241</td><td>241</td><td>-0.00848114</td><td>-0.009318300</td><td>-0.01264550</td><td>-0.024159000</td><td> 0.004272330</td><td>-0.00742088</td><td>-0.021371400</td><td> 0.006188810</td><td>⋯</td><td>-0.020598800</td><td>-0.00816347</td><td>-0.00669854</td><td>-0.002650730</td><td>-0.000480313</td><td> 0.022766400</td><td>-0.011031400</td><td> 0.02638610</td><td> 0.041146100</td><td>-0.05040010</td></tr>
	<tr><td>251</td><td>251</td><td>-0.00321513</td><td>-0.001022570</td><td>-0.01629190</td><td>-0.053413500</td><td>-0.008869390</td><td>-0.02919450</td><td>-0.041730300</td><td>-0.052242800</td><td>⋯</td><td> 0.060614500</td><td> 0.02645980</td><td> 0.03396130</td><td>-0.017617200</td><td>-0.018829300</td><td>-0.046824100</td><td>-0.017043900</td><td>-0.06253150</td><td> 0.051962100</td><td> 0.00823284</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>10936</td><td>10936</td><td> 0.192199000</td><td>-1.21779e-01</td><td> 0.014781500</td><td> 2.38885e-02</td><td>-0.01441160</td><td>-0.041398000</td><td> 0.023953400</td><td> 0.02070130</td><td>⋯</td><td>-0.009883150</td><td> 0.00322184</td><td>-0.006657390</td><td>-0.05924230</td><td>-0.022464900</td><td>-0.21990500</td><td> 0.05018110</td><td> 0.02376400</td><td>-0.128982000</td><td>-0.18018800</td></tr>
	<tr><td>10962</td><td>10962</td><td>-0.004675960</td><td>-7.31984e-03</td><td>-0.002711820</td><td>-1.08159e-02</td><td>-0.00516487</td><td>-0.000516792</td><td> 0.005821840</td><td>-0.00312742</td><td>⋯</td><td>-0.053152600</td><td>-0.01283800</td><td>-0.030819800</td><td> 0.01503450</td><td> 0.018971300</td><td> 0.01810600</td><td> 0.05058500</td><td> 0.02571450</td><td>-0.012599500</td><td> 0.00730739</td></tr>
	<tr><td>10967</td><td>10967</td><td>-0.011491700</td><td>-4.59450e-03</td><td>-0.017790500</td><td>-1.36199e-03</td><td>-0.00418879</td><td>-0.020403500</td><td> 0.005149450</td><td> 0.00501159</td><td>⋯</td><td>-0.054938500</td><td> 0.00563654</td><td> 0.006627810</td><td> 0.05799610</td><td>-0.006413230</td><td> 0.00368716</td><td> 0.02052390</td><td> 0.00139196</td><td>-0.013971400</td><td> 0.00728833</td></tr>
	<tr><td>10977</td><td>10977</td><td>-0.011301200</td><td>-2.06668e-03</td><td> 0.002961410</td><td>-2.61940e-03</td><td> 0.01182870</td><td> 0.004583340</td><td>-0.021214300</td><td> 0.02589610</td><td>⋯</td><td>-0.042936300</td><td>-0.01176170</td><td> 0.030962000</td><td>-0.02672850</td><td> 0.013566000</td><td>-0.00777501</td><td>-0.00175320</td><td> 0.03284810</td><td> 0.058455900</td><td>-0.02569390</td></tr>
	<tr><td>10984</td><td>10984</td><td>-0.011978600</td><td>-3.23065e-03</td><td> 0.007256620</td><td> 1.32392e-02</td><td>-0.00783351</td><td> 0.051677100</td><td> 0.002973850</td><td>-0.01786740</td><td>⋯</td><td> 0.029385200</td><td> 0.00510882</td><td>-0.024039200</td><td>-0.01874260</td><td>-0.008857520</td><td>-0.01992440</td><td> 0.01922710</td><td> 0.02887310</td><td> 0.003494130</td><td>-0.00602685</td></tr>
	<tr><td>11015</td><td>11015</td><td>-0.006894840</td><td>-6.58623e-03</td><td> 0.011427900</td><td> 1.94331e-02</td><td> 0.00418976</td><td> 0.032980900</td><td>-0.035573400</td><td> 0.01164990</td><td>⋯</td><td> 0.055893000</td><td>-0.00249302</td><td> 0.046919000</td><td>-0.05949310</td><td>-0.000597924</td><td>-0.00840922</td><td> 0.02007410</td><td> 0.04578450</td><td>-0.028155400</td><td>-0.00925084</td></tr>
	<tr><td>11053</td><td>11053</td><td>-0.008526660</td><td>-4.60711e-03</td><td> 0.009316500</td><td> 1.97597e-02</td><td> 0.00454710</td><td>-0.026175200</td><td>-0.044525300</td><td> 0.04552800</td><td>⋯</td><td> 0.046319700</td><td> 0.03918930</td><td> 0.044371700</td><td> 0.06591800</td><td> 0.000383875</td><td> 0.01143710</td><td> 0.03452430</td><td>-0.03654790</td><td>-0.005471120</td><td>-0.04007420</td></tr>
	<tr><td>11069</td><td>11069</td><td>-0.008577040</td><td> 3.55202e-03</td><td>-0.001171780</td><td>-4.75477e-06</td><td>-0.00936589</td><td>-0.038486200</td><td>-0.020122500</td><td>-0.06094120</td><td>⋯</td><td> 0.025537900</td><td>-0.02199680</td><td> 0.074407000</td><td>-0.00417962</td><td> 0.027125400</td><td>-0.02607070</td><td> 0.02730820</td><td>-0.00579564</td><td> 0.035290800</td><td> 0.00888144</td></tr>
	<tr><td>11076</td><td>11076</td><td>-0.009207010</td><td>-3.29066e-03</td><td> 0.008499500</td><td> 3.67877e-03</td><td> 0.00478666</td><td>-0.009053030</td><td>-0.003448130</td><td> 0.06725430</td><td>⋯</td><td>-0.032421900</td><td>-0.03549710</td><td> 0.007698970</td><td>-0.06273900</td><td> 0.007267200</td><td> 0.03091170</td><td> 0.02713880</td><td> 0.08387770</td><td> 0.001821350</td><td>-0.00817811</td></tr>
	<tr><td>11135</td><td>11135</td><td> 0.020058400</td><td> 4.30007e-02</td><td> 0.022473500</td><td> 7.38412e-03</td><td>-0.02401760</td><td>-0.031742600</td><td>-0.029018600</td><td> 0.03229260</td><td>⋯</td><td> 0.017091600</td><td>-0.04799250</td><td>-0.001228870</td><td>-0.06561860</td><td> 0.034261200</td><td> 0.05038870</td><td>-0.06292640</td><td>-0.04337150</td><td>-0.033266000</td><td>-0.06188330</td></tr>
	<tr><td>11388</td><td>11388</td><td>-0.011078400</td><td>-2.77155e-03</td><td> 0.007895890</td><td>-4.59851e-03</td><td>-0.00799603</td><td>-0.016429000</td><td> 0.013055300</td><td> 0.03212730</td><td>⋯</td><td> 0.008152360</td><td> 0.02449630</td><td> 0.007470700</td><td>-0.00482444</td><td>-0.002022360</td><td> 0.00708970</td><td> 0.00523126</td><td> 0.01370950</td><td> 0.038273200</td><td> 0.00120693</td></tr>
	<tr><td>11396</td><td>11396</td><td>-0.004497760</td><td>-8.60218e-04</td><td>-0.011930500</td><td>-5.46915e-02</td><td>-0.00792371</td><td>-0.022179600</td><td>-0.048340400</td><td> 0.05409370</td><td>⋯</td><td>-0.047029900</td><td>-0.01304940</td><td>-0.034625900</td><td>-0.01206940</td><td>-0.016569200</td><td> 0.02329750</td><td> 0.02238070</td><td>-0.03960180</td><td> 0.034664700</td><td> 0.02139500</td></tr>
	<tr><td>11397</td><td>11397</td><td>-0.013367100</td><td>-1.75078e-03</td><td> 0.019612700</td><td> 1.94201e-02</td><td> 0.00160817</td><td>-0.033233200</td><td> 0.051239600</td><td> 0.03394310</td><td>⋯</td><td>-0.028413500</td><td> 0.01451720</td><td>-0.008638270</td><td>-0.04327210</td><td> 0.007880220</td><td> 0.02066400</td><td> 0.00709681</td><td> 0.02662350</td><td>-0.002449890</td><td>-0.03784340</td></tr>
	<tr><td>11425</td><td>11425</td><td>-0.007650230</td><td>-6.93880e-03</td><td> 0.010349200</td><td> 1.23778e-02</td><td> 0.00384237</td><td>-0.023796400</td><td>-0.007666500</td><td> 0.00939161</td><td>⋯</td><td> 0.079037400</td><td>-0.00734664</td><td>-0.018063000</td><td> 0.05781490</td><td>-0.010909600</td><td> 0.01995370</td><td>-0.00999942</td><td>-0.01188510</td><td>-0.003837410</td><td> 0.06226800</td></tr>
	<tr><td>11441</td><td>11441</td><td>-0.009893980</td><td>-2.33273e-03</td><td> 0.002324280</td><td>-4.27340e-03</td><td>-0.01592630</td><td>-0.025024800</td><td>-0.010839500</td><td> 0.03767540</td><td>⋯</td><td>-0.021090200</td><td> 0.03203370</td><td>-0.006628100</td><td>-0.02706470</td><td>-0.000570901</td><td> 0.02155580</td><td> 0.01338160</td><td> 0.01483090</td><td>-0.029400100</td><td>-0.04796310</td></tr>
	<tr><td>11458</td><td>11458</td><td>-0.012870900</td><td>-4.63183e-03</td><td> 0.000472721</td><td> 8.27262e-03</td><td> 0.00177484</td><td>-0.024239400</td><td>-0.006930730</td><td>-0.03505970</td><td>⋯</td><td> 0.021579700</td><td>-0.05769070</td><td>-0.000757448</td><td>-0.01854510</td><td> 0.042715800</td><td>-0.01023350</td><td> 0.02850570</td><td> 0.01098940</td><td>-0.002778310</td><td> 0.00937549</td></tr>
	<tr><td>11493</td><td>11493</td><td>-0.007011380</td><td>-1.84579e-04</td><td> 0.006584320</td><td> 1.66167e-02</td><td> 0.01382400</td><td>-0.022373100</td><td>-0.022429100</td><td> 0.05440250</td><td>⋯</td><td> 0.120426000</td><td> 0.03180810</td><td> 0.037522300</td><td> 0.02083650</td><td>-0.062404600</td><td>-0.02686040</td><td>-0.04938200</td><td>-0.01505490</td><td>-0.000790275</td><td>-0.00403634</td></tr>
	<tr><td>11519</td><td>11519</td><td> 0.002960570</td><td> 2.68257e-03</td><td>-0.029608900</td><td>-7.46827e-02</td><td> 0.01757660</td><td> 0.004000530</td><td> 0.000158212</td><td> 0.00477255</td><td>⋯</td><td>-0.043145600</td><td>-0.06457940</td><td>-0.029691700</td><td>-0.00997568</td><td> 0.000040224</td><td> 0.00909673</td><td>-0.00293489</td><td>-0.00934990</td><td>-0.077850600</td><td> 0.04556860</td></tr>
	<tr><td>11525</td><td>11525</td><td>-0.000899349</td><td> 1.79833e-02</td><td> 0.031515000</td><td> 4.66881e-02</td><td> 0.01313700</td><td>-0.045685000</td><td>-0.084504700</td><td>-0.08657600</td><td>⋯</td><td>-0.061848400</td><td>-0.06218010</td><td>-0.039155900</td><td> 0.00169151</td><td> 0.033602700</td><td> 0.05550260</td><td> 0.07227110</td><td> 0.00595697</td><td>-0.004855320</td><td>-0.00183803</td></tr>
	<tr><td>11526</td><td>11526</td><td>-0.012813700</td><td>-7.98356e-03</td><td>-0.000921799</td><td> 2.26999e-02</td><td> 0.00231765</td><td>-0.008596100</td><td> 0.012463900</td><td> 0.01913980</td><td>⋯</td><td> 0.006848640</td><td> 0.00580370</td><td>-0.006117520</td><td>-0.01926040</td><td>-0.008190520</td><td> 0.04534460</td><td>-0.03446500</td><td>-0.01815050</td><td> 0.013650300</td><td>-0.01362000</td></tr>
	<tr><td>11528</td><td>11528</td><td>-0.012371900</td><td>-3.72862e-03</td><td> 0.002948210</td><td> 9.63867e-03</td><td>-0.00678087</td><td>-0.011237300</td><td> 0.023412700</td><td> 0.00434507</td><td>⋯</td><td>-0.037989000</td><td> 0.01450220</td><td> 0.067127500</td><td> 0.00672536</td><td>-0.009651970</td><td>-0.02913130</td><td> 0.00210938</td><td>-0.01930910</td><td>-0.000139403</td><td> 0.01341480</td></tr>
	<tr><td>11531</td><td>11531</td><td>-0.011288400</td><td>-5.94892e-03</td><td> 0.006610880</td><td> 1.72562e-02</td><td> 0.00866793</td><td>-0.039839900</td><td>-0.018248200</td><td> 0.01334560</td><td>⋯</td><td>-0.036971400</td><td>-0.01274840</td><td> 0.010981300</td><td>-0.03858470</td><td> 0.024759100</td><td>-0.00106921</td><td>-0.03010670</td><td>-0.00144590</td><td>-0.013109100</td><td> 0.01531040</td></tr>
	<tr><td>11562</td><td>11562</td><td>-0.008730700</td><td>-7.99532e-03</td><td> 0.006836860</td><td>-5.09707e-03</td><td> 0.00320566</td><td>-0.032107300</td><td>-0.027381800</td><td> 0.04714260</td><td>⋯</td><td> 0.014993500</td><td> 0.00259318</td><td>-0.006275940</td><td>-0.05706420</td><td>-0.001731370</td><td> 0.00856210</td><td> 0.01404350</td><td>-0.03031540</td><td>-0.010866200</td><td>-0.00393769</td></tr>
	<tr><td>11584</td><td>11584</td><td>-0.006362800</td><td> 1.05591e-03</td><td> 0.010858900</td><td> 2.44756e-02</td><td>-0.01701260</td><td> 0.061087200</td><td> 0.010992100</td><td> 0.04121230</td><td>⋯</td><td> 0.047344700</td><td>-0.03928830</td><td> 0.041021500</td><td> 0.05189450</td><td>-0.009976000</td><td> 0.00249817</td><td> 0.00290588</td><td> 0.01662550</td><td> 0.014133000</td><td>-0.02094810</td></tr>
	<tr><td>11641</td><td>11641</td><td>-0.009977770</td><td>-5.79947e-03</td><td>-0.010518500</td><td> 6.89173e-03</td><td> 0.01300410</td><td>-0.022218300</td><td>-0.011946300</td><td>-0.08912890</td><td>⋯</td><td> 0.029790000</td><td> 0.02744520</td><td>-0.048400900</td><td>-0.01832100</td><td> 0.004571160</td><td> 0.03408650</td><td>-0.02956850</td><td>-0.02508840</td><td> 0.053773000</td><td>-0.02705210</td></tr>
	<tr><td>11698</td><td>11698</td><td>-0.009587820</td><td>-4.41876e-03</td><td>-0.009336360</td><td>-3.29014e-03</td><td>-0.00147518</td><td>-0.019141200</td><td> 0.006812520</td><td>-0.04861050</td><td>⋯</td><td> 0.039945700</td><td> 0.00769957</td><td> 0.012009700</td><td> 0.01303730</td><td> 0.027260600</td><td>-0.01942780</td><td> 0.01997140</td><td> 0.03009660</td><td> 0.004389110</td><td>-0.02717940</td></tr>
	<tr><td>11753</td><td>11753</td><td>-0.013639200</td><td>-1.17043e-02</td><td>-0.000855109</td><td> 1.69574e-02</td><td>-0.01664600</td><td> 0.106480000</td><td>-0.016796100</td><td>-0.00977177</td><td>⋯</td><td> 0.000839177</td><td> 0.04744920</td><td>-0.023221000</td><td> 0.01150480</td><td>-0.034613700</td><td>-0.01948790</td><td> 0.01182420</td><td>-0.01044920</td><td>-0.005756000</td><td>-0.04476280</td></tr>
	<tr><td>11795</td><td>11795</td><td>-0.003828500</td><td>-6.68915e-05</td><td> 0.016649400</td><td> 2.25705e-02</td><td> 0.00590070</td><td>-0.025347700</td><td>-0.016259700</td><td> 0.02598540</td><td>⋯</td><td>-0.035563200</td><td>-0.03308420</td><td>-0.006697990</td><td> 0.02139710</td><td>-0.004562050</td><td> 0.01179730</td><td> 0.05813120</td><td>-0.01254420</td><td> 0.042140500</td><td> 0.01582580</td></tr>
	<tr><td>11802</td><td>11802</td><td>-0.003638700</td><td> 1.86601e-02</td><td> 0.017731200</td><td> 3.59520e-02</td><td> 0.00745392</td><td>-0.051820500</td><td>-0.069269500</td><td>-0.06772790</td><td>⋯</td><td>-0.030866000</td><td>-0.02864320</td><td> 0.001669760</td><td>-0.01465850</td><td>-0.002653950</td><td> 0.04898730</td><td> 0.04365200</td><td> 0.00751692</td><td>-0.017578800</td><td> 0.02448530</td></tr>
	<tr><td>11819</td><td>11819</td><td>-0.012335800</td><td>-9.89583e-05</td><td> 0.022142500</td><td> 1.29498e-02</td><td>-0.00198189</td><td> 0.019274200</td><td>-0.017522300</td><td> 0.03581650</td><td>⋯</td><td>-0.037033200</td><td> 0.03027750</td><td> 0.026589900</td><td>-0.03930310</td><td> 0.024352700</td><td>-0.00806471</td><td>-0.01701110</td><td> 0.00440610</td><td>-0.048689600</td><td>-0.00149517</td></tr>
</tbody>
</table>




```R
# Read eigenvalues file
eigenval <- scan("merged_allchips_common_pca.eigenval")

# Calculate proportion variance explained per PC
var_exp <- eigenval / sum(eigenval) * 100

# Check first few variances
head(var_exp)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>28.7817568177969</li><li>12.5810999241619</li><li>5.09010343330398</li><li>4.36439806498501</li><li>4.0240768453216</li><li>3.57761213888523</li></ol>




```R
library(ggplot2)

ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6) +
  xlab(paste0("PC1 (", round(var_exp[1], 2), "% variance)")) +
  ylab(paste0("PC2 (", round(var_exp[2], 2), "% variance)")) +
  theme_minimal() +
  ggtitle("PCA Plot: PC1 vs PC2") +
  theme(plot.title = element_text(hjust = 0.5))

```


    
![png](output_62_0.png)
    



```R
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata$chip[metadata$chip == ""] <- "Unknown"
head(metadata)
colnames(metadata)
```


<table class="dataframe">
<caption>A data.frame: 6 × 6</caption>
<thead>
	<tr><th></th><th scope=col>user</th><th scope=col>build</th><th scope=col>chip</th><th scope=col>chip_version</th><th scope=col>inferred_sex</th><th scope=col>source</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>10645</td><td>37</td><td>Unknown       </td><td>  </td><td>Female</td><td>AncestryDNA             </td></tr>
	<tr><th scope=row>2</th><td>10542</td><td>37</td><td>HTS iSelect HD</td><td>v4</td><td>Female</td><td>23andMe                 </td></tr>
	<tr><th scope=row>3</th><td>10653</td><td>37</td><td>Unknown       </td><td>  </td><td>Female</td><td>AncestryDNA, AncestryDNA</td></tr>
	<tr><th scope=row>4</th><td>10106</td><td>37</td><td>Unknown       </td><td>  </td><td>Female</td><td>AncestryDNA             </td></tr>
	<tr><th scope=row>5</th><td>10559</td><td>37</td><td>Illumina GSAs </td><td>v5</td><td>Male  </td><td>23andMe                 </td></tr>
	<tr><th scope=row>6</th><td>10828</td><td>37</td><td>Illumina GSAs </td><td>v5</td><td>Female</td><td>23andMe                 </td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'user'</li><li>'build'</li><li>'chip'</li><li>'chip_version'</li><li>'inferred_sex'</li><li>'source'</li></ol>




```R
colnames(metadata) <- c("FID", "build", "chip", "chip_version", "inferred_sex", "source")
```


```R
pcs_meta <- merge(pcs, metadata, by = c("FID"))
head(pcs_meta)
```


<table class="dataframe">
<caption>A data.frame: 6 × 27</caption>
<thead>
	<tr><th></th><th scope=col>FID</th><th scope=col>IID</th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>⋯</th><th scope=col>PC16</th><th scope=col>PC17</th><th scope=col>PC18</th><th scope=col>PC19</th><th scope=col>PC20</th><th scope=col>build</th><th scope=col>chip</th><th scope=col>chip_version</th><th scope=col>inferred_sex</th><th scope=col>source</th></tr>
	<tr><th></th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td> 8</td><td> 8</td><td>-0.00695248</td><td>-0.000978118</td><td>-0.01889710</td><td>-0.04086600</td><td>-0.00226661</td><td> 0.0351734</td><td>-0.000920792</td><td>-0.0516592</td><td>⋯</td><td> 0.02755300</td><td>-0.01868250</td><td> 0.01440690</td><td>-0.00737957</td><td> 0.02512950</td><td>37</td><td>OmniExpress plus</td><td>v3</td><td>Female</td><td>23andMe     </td></tr>
	<tr><th scope=row>2</th><td>11</td><td>11</td><td> 0.00638909</td><td>-0.000885972</td><td>-0.06231530</td><td>-0.11414400</td><td>-0.02321030</td><td> 0.0664040</td><td>-0.011471100</td><td>-0.0501895</td><td>⋯</td><td> 0.00625582</td><td>-0.01436440</td><td> 0.01456890</td><td> 0.01708750</td><td> 0.02230760</td><td>37</td><td>OmniExpress plus</td><td>v3</td><td>Male  </td><td>23andMe     </td></tr>
	<tr><th scope=row>3</th><td>14</td><td>14</td><td>-0.01231010</td><td>-0.004080500</td><td> 0.00722458</td><td> 0.00892257</td><td> 0.00432515</td><td>-0.0216232</td><td> 0.027711000</td><td> 0.0421512</td><td>⋯</td><td> 0.00667454</td><td>-0.01709200</td><td>-0.02619020</td><td>-0.01470450</td><td>-0.01178390</td><td>37</td><td>Unknown         </td><td>  </td><td>Female</td><td>23andMe     </td></tr>
	<tr><th scope=row>4</th><td>16</td><td>16</td><td>-0.01386140</td><td>-0.002910680</td><td>-0.00355122</td><td> 0.01436510</td><td>-0.01112460</td><td> 0.0391532</td><td>-0.002504390</td><td>-0.0215492</td><td>⋯</td><td>-0.07659920</td><td>-0.01784900</td><td>-0.00166217</td><td> 0.00407523</td><td> 0.02687590</td><td>37</td><td>HTS iSelect HD  </td><td>v4</td><td>Male  </td><td>23andMe     </td></tr>
	<tr><th scope=row>5</th><td>17</td><td>17</td><td>-0.01186810</td><td>-0.001851780</td><td> 0.01857800</td><td> 0.02757830</td><td>-0.00194662</td><td>-0.0217239</td><td>-0.013766200</td><td> 0.0398188</td><td>⋯</td><td>-0.03242860</td><td> 0.04441440</td><td> 0.00827015</td><td> 0.03226100</td><td>-0.01576650</td><td>37</td><td>Unknown         </td><td>  </td><td>Female</td><td>23andMe     </td></tr>
	<tr><th scope=row>6</th><td>26</td><td>26</td><td>-0.00637205</td><td> 0.002338360</td><td> 0.01158400</td><td> 0.01872460</td><td> 0.01837330</td><td>-0.0293547</td><td>-0.005018150</td><td>-0.0169588</td><td>⋯</td><td> 0.03348180</td><td> 0.00679015</td><td>-0.02563430</td><td> 0.02159630</td><td> 0.00792364</td><td>37</td><td>Unknown         </td><td>  </td><td>Male  </td><td>vcf, 23andMe</td></tr>
</tbody>
</table>




```R
library(ggplot2)

ggplot(pcs_meta, aes(x = PC1, y = PC2, color = chip)) +
  geom_point(alpha = 0.7, size = 2) +
  xlab(paste0("PC1 (", round(var_exp[1], 2), "% variance)")) +
  ylab(paste0("PC2 (", round(var_exp[2], 2), "% variance)")) +
  theme_minimal() +
  ggtitle("PCA Plot Colored by Chip") +
  theme(plot.title = element_text(hjust = 0.5))

```


    
![png](output_66_0.png)
    



```R
ggplot(pcs_meta, aes(x = PC1, y = PC2, color = inferred_sex)) +
  geom_point(alpha = 0.7, size = 2) +
  xlab(paste0("PC1 (", round(var_exp[1], 2), "% variance)")) +
  ylab(paste0("PC2 (", round(var_exp[2], 2), "% variance)")) +
  theme_minimal() +
  ggtitle("PCA Plot Colored by Sex") +
  theme(plot.title = element_text(hjust = 0.5))

```


    
![png](output_67_0.png)
    



```R
phenotype <- read.table("height.txt", header = TRUE)
colnames(phenotype) <- c("FID", "Height")
pcs_meta_pheno <- merge(pcs_meta, phenotype, by = c("FID"))
head(pcs_meta_pheno)
```


<table class="dataframe">
<caption>A data.frame: 6 × 28</caption>
<thead>
	<tr><th></th><th scope=col>FID</th><th scope=col>IID</th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>⋯</th><th scope=col>PC17</th><th scope=col>PC18</th><th scope=col>PC19</th><th scope=col>PC20</th><th scope=col>build</th><th scope=col>chip</th><th scope=col>chip_version</th><th scope=col>inferred_sex</th><th scope=col>source</th><th scope=col>Height</th></tr>
	<tr><th></th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td> 8</td><td> 8</td><td>-0.00695248</td><td>-0.000978118</td><td>-0.01889710</td><td>-0.04086600</td><td>-0.00226661</td><td> 0.0351734</td><td>-0.000920792</td><td>-0.0516592</td><td>⋯</td><td>-0.01868250</td><td> 0.01440690</td><td>-0.00737957</td><td> 0.02512950</td><td>37</td><td>OmniExpress plus</td><td>v3</td><td>Female</td><td>23andMe     </td><td>171.00</td></tr>
	<tr><th scope=row>2</th><td>11</td><td>11</td><td> 0.00638909</td><td>-0.000885972</td><td>-0.06231530</td><td>-0.11414400</td><td>-0.02321030</td><td> 0.0664040</td><td>-0.011471100</td><td>-0.0501895</td><td>⋯</td><td>-0.01436440</td><td> 0.01456890</td><td> 0.01708750</td><td> 0.02230760</td><td>37</td><td>OmniExpress plus</td><td>v3</td><td>Male  </td><td>23andMe     </td><td>177.80</td></tr>
	<tr><th scope=row>3</th><td>14</td><td>14</td><td>-0.01231010</td><td>-0.004080500</td><td> 0.00722458</td><td> 0.00892257</td><td> 0.00432515</td><td>-0.0216232</td><td> 0.027711000</td><td> 0.0421512</td><td>⋯</td><td>-0.01709200</td><td>-0.02619020</td><td>-0.01470450</td><td>-0.01178390</td><td>37</td><td>Unknown         </td><td>  </td><td>Female</td><td>23andMe     </td><td>155.00</td></tr>
	<tr><th scope=row>4</th><td>16</td><td>16</td><td>-0.01386140</td><td>-0.002910680</td><td>-0.00355122</td><td> 0.01436510</td><td>-0.01112460</td><td> 0.0391532</td><td>-0.002504390</td><td>-0.0215492</td><td>⋯</td><td>-0.01784900</td><td>-0.00166217</td><td> 0.00407523</td><td> 0.02687590</td><td>37</td><td>HTS iSelect HD  </td><td>v4</td><td>Male  </td><td>23andMe     </td><td>183.00</td></tr>
	<tr><th scope=row>5</th><td>17</td><td>17</td><td>-0.01186810</td><td>-0.001851780</td><td> 0.01857800</td><td> 0.02757830</td><td>-0.00194662</td><td>-0.0217239</td><td>-0.013766200</td><td> 0.0398188</td><td>⋯</td><td> 0.04441440</td><td> 0.00827015</td><td> 0.03226100</td><td>-0.01576650</td><td>37</td><td>Unknown         </td><td>  </td><td>Female</td><td>23andMe     </td><td>172.72</td></tr>
	<tr><th scope=row>6</th><td>26</td><td>26</td><td>-0.00637205</td><td> 0.002338360</td><td> 0.01158400</td><td> 0.01872460</td><td> 0.01837330</td><td>-0.0293547</td><td>-0.005018150</td><td>-0.0169588</td><td>⋯</td><td> 0.00679015</td><td>-0.02563430</td><td> 0.02159630</td><td> 0.00792364</td><td>37</td><td>Unknown         </td><td>  </td><td>Male  </td><td>vcf, 23andMe</td><td>188.00</td></tr>
</tbody>
</table>




```R
pcs_meta_pheno$Height_bin <- cut(
  pcs_meta_pheno$Height,
  breaks = seq(140, 200, by = 10),  # Adjust as needed
  right = FALSE,
  labels = paste0(seq(140, 190, by = 10), "-", seq(149, 199, by = 10))
)
```


```R
ggplot(pcs_meta_pheno, aes(x = PC1, y = PC2, color = Height_bin)) +
  geom_point(alpha = 0.7, size = 2) +
  xlab(paste0("PC1 (", round(var_exp[1], 2), "% variance)")) +
  ylab(paste0("PC2 (", round(var_exp[2], 2), "% variance)")) +
  theme_minimal() +
  ggtitle("PCA Plot Colored by Height") +
  theme(plot.title = element_text(hjust = 0.5))
```


    
![png](output_70_0.png)
    


## File preparation for association analysis is done in another file "last chance". Disregard all of the next lines


```R

```


```R
# Read fam file
fam <- read.table("merged_allchips_common.fam", header=FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# Extract FID and IID from fam
ids <- fam[, c("FID", "IID")]

# Read height.txt (two columns: IID and HEIGHT)
height <- read.table("height.txt", header=FALSE)
colnames(height) <- c("IID", "HEIGHT")

# Since FID == IID in fam, create FID in height for merge
height$FID <- height$IID

# Reorder columns to match for merging: FID, IID, HEIGHT
height <- height[, c("FID", "IID", "HEIGHT")]

# Merge on FID and IID
pheno_merged <- merge(ids, height, by=c("FID", "IID"), all.x=TRUE)

# Check for missing phenotypes
missing_count <- sum(is.na(pheno_merged$HEIGHT))
if(missing_count > 0){
  cat("Warning: There are", missing_count, "samples missing phenotype data!\n")
}

# Rename HEIGHT to PHENO for PLINK
colnames(pheno_merged)[3] <- "PHENO"

# Write phenotype.txt
write.table(pheno_merged, file="phenotype.txt", quote=FALSE, row.names=FALSE)

```


```R
# Read ids.txt (contains FID and IID)
ids <- read.table("ids.txt", header=FALSE, col.names=c("FID", "IID"))

# Read your full phenotype.txt (with FID, IID, PHENO)
phenotype <- read.table("phenotype.txt", header=TRUE)

# Filter phenotype to keep only rows present in ids.txt
phenotype_filtered <- merge(ids, phenotype, by=c("FID", "IID"))

# Write filtered phenotype file (no row names, no quotes)
write.table(phenotype_filtered, file="phenotype_filtered.txt", 
            quote=FALSE, row.names=FALSE)

```


```R
# Read phenotype_filtered.txt
phenotype_filtered <- read.table("phenotype_filtered.txt", header=TRUE)

# Check duplicates by FID and IID
duplicates <- phenotype_filtered[duplicated(phenotype_filtered[, c("FID", "IID")]) |
                                 duplicated(phenotype_filtered[, c("FID", "IID")], fromLast=TRUE), ]

# Show duplicates if any
if(nrow(duplicates) > 0){
  print("Duplicate FID-IID pairs found:")
  print(duplicates)
} else {
  print("No duplicates found in phenotype_filtered.txt")
}
```

    [1] "No duplicates found in phenotype_filtered.txt"


covariates file


```R
# Read fam file (has SEX info)
fam <- read.table("merged_allchips_common.fam", header=FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# Read ids.txt (filtered sample list)
ids <- read.table("ids.txt", header=FALSE, col.names=c("FID", "IID"))

# Merge to keep only covariates for IDs in ids.txt
covariates <- merge(ids, fam[, c("FID", "IID", "SEX")], by=c("FID", "IID"))

# Optional: check if any samples are missing SEX info (usually 1 or 2)
missing_sex <- sum(is.na(covariates$SEX))
if(missing_sex > 0) {
  warning(paste(missing_sex, "samples have missing SEX info"))
}
missing_sex
```


0



```R
# Write covariates.txt
write.table(covariates, "covariates.txt", quote=FALSE, row.names=FALSE)
```


```R

```


```R
# Check for NA or invalid p-values
summary(gwas$P)
anyNA(gwas$P)         # TRUE if there are NAs
any(!is.finite(gwas$P)) # TRUE if Inf or -Inf exist

# Remove rows with missing or invalid p-values before plotting
gwas_clean <- gwas[!is.na(gwas$P) & is.finite(gwas$P), ]

# Then try QQ plot again:
library(qqman)
qq(gwas_clean$P, main="QQ plot of GWAS P-values")

# And Manhattan plot:
manhattan(gwas_clean, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot", genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5))

```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0.000   0.000   0.000   0.233   0.454   1.000    4458 



TRUE



TRUE



    
![png](output_80_3.png)
    



    
![png](output_80_4.png)
    



```R
summary(results_clean$P)
table(is.finite(results_clean$P))

```


        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.000000 0.003382 0.005890 0.132473 0.154500 1.000000 



    
      TRUE 
    261835 



```R
summary(results$P)
table(is.finite(results$P))
```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0.000   0.003   0.006   0.132   0.154   1.000    8916 



    
     FALSE   TRUE 
      8916 261835 



```R
library(qqman)
qq(results_clean$P)

```


    
![png](output_83_0.png)
    



```R

```


```R

```


```R

```


```R
library(ggplot2)

ggplot(all_lmiss, aes(x = F_MISS, fill = chip)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
  theme_minimal() +
  labs(title = "SNP Missingness per Chip Type",
       x = "Fraction Missing (F_MISS)",
       y = "Number of SNPs") +
  theme(legend.position = "right")
```
