library("dplyr")
library("tidyr")
X <- read.table("data/myc_zinc_splitroot_robbe.tsv",dec=",",sep="\t",header=T)

# Small check:
table(X$batch,X$condition)

# First, make a simplified version of the dataframe, for looking at masses
# Subset dataframe with first occurrence of each unique sample-strain combination
# We take the means of the mycorrhization percentages, and count the number of "y"

x <- X %>%
  group_by(sample, strain) %>%
  summarise(
    condition = first(condition),  # Retain condition for each sample
    batch = first(batch),  # Retain batch for each sample
    root_mass = first(root_mass),  # Keeping first occurrence (assuming same across triplicates)
    shoot_mass = first(shoot_mass),  # Same assumption
    mycorrhization_percentage = mean(mycorrhization_percentage, na.rm = TRUE),  # Mean mycorrhization
    root_hairs_count = sum(root_hairs == "y", na.rm = TRUE),  # Count occurrences of "y"
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = strain,
    values_from = c(root_mass, mycorrhization_percentage, root_hairs_count),
    names_glue = "{.value}_{strain}"
  ) %>%
  mutate(
    across(c(condition, batch, sample), as.factor),  # Ensure categorical variables are factors
    across(where(is.numeric), as.numeric)  # Ensure numeric columns remain numeric
  ) %>%
  as.data.frame()

# Next, B56 and B51 turn out to be an outlier on every level. They barely had roots. So we just kick them out. 
x <- x[x$sample!="B56",]
x <- x[x$sample!="B51",]

# relevel x$condition for a more logical order:
x$condition <- factor(x$condition, levels = c("control","50µM","200µM"))

x$root_logratio           <- log2(x$root_mass_LS4 / x$root_mass_HH1)
x$nonmycorrhized_logratio <- log2((100-x$mycorrhization_percentage_LS4)/(100-x$mycorrhization_percentage_HH1))
x$mycorrhization_diff     <- x$mycorrhization_percentage_LS4 - x$mycorrhization_percentage_HH1
x$root_proportion_LS4     <- x$root_mass_LS4 / (x$root_mass_LS4 + x$root_mass_HH1 + x$shoot_mass)
x$root_proportion_HH1     <- x$root_mass_HH1 / (x$root_mass_LS4 + x$root_mass_HH1 + x$shoot_mass)
x$shoot_proportion        <- x$shoot_mass    / (x$shoot_mass + x$root_mass_HH1 + x$root_mass_LS4)
x$condition_numeric       <- c(0,50,200)[x$condition]

# Interesting observations:

par(mfrow=c(1,2))
plot(x$mycorrhization_percentage_HH1 ~ x$condition,ylab="% mycorrhiza",xlab="condition",main="HH1 (sensitive strain)")
plot(x$mycorrhization_percentage_LS4 ~ x$condition,ylab="% mycorrhiza",xlab="condition",main="LS4 (tolerant strain)")

lm1a <- lm(x$mycorrhization_percentage_LS4 ~ x$condition + x$mycorrhization_percentage_HH1)
lm1b <- lm(x$mycorrhization_percentage_LS4 ~ x$condition_binary + x$mycorrhization_percentage_HH1)
lm2a <- lm(x$mycorrhization_percentage_HH1 ~ x$condition + x$mycorrhization_percentage_LS4)
lm2b <- lm(x$mycorrhization_percentage_HH1 ~ x$condition_binary + x$mycorrhization_percentage_LS4)

anova(lm1a) # etc
# I also tested interactions but that went nowhere

# We should check if the outlier is not something strange
plot(x$shoot_mass ~ x$condition)
plot(x$shoot_mass ~ x$condition_binary)




# First very broadly looking at differences between strains:
boxplot(x$root_mass~x$strain, main = "root mass")
boxplot(x$root_proportion~x$strain, main = "root proportion")
boxplot(x$shoot_proportion~x$strain, main = "shoot proportion")

# There clearly is an outlier that does not seem to have roots. 


train_index <- createDataPartition(x_rf$condition, p = 0.7, list = FALSE)
train_data <- x_rf[train_index, ]
test_data  <- x_rf[-train_index, ]

