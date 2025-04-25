library("dplyr")
library("tidyr")
library("caret")
library("randomForest")
library("piecewiseSEM")
library("lme4")        
library("multcompView")
library("lavaan")
library("semPlot")


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
x$condition_binary        <- as.factor(c("control","test")[as.factor(x$condition!="control")])

# Interesting observations:
par(mfrow=c(1,2))

# Batch effects?
boxplot(x$mycorrhization_percentage_LS4 ~ x$batch,col=c("green","orange","red")[x$condition[!duplicated(x$batch)]],main="LS4",xlab="batch",ylab="myc%",ylim=c(0,100))
boxplot(x$mycorrhization_percentage_HH1 ~ x$batch,col=c("green","orange","red")[x$condition[!duplicated(x$batch)]],main="HH1",xlab="batch",ylab="myc%",ylim=c(0,100))
legend("topright",legend=levels(x$condition),col=c("green","orange","red"),pch=15)



# The percentage of mycorrhiza seems to be stable in the sensitive strain, while being higher in the LS4 strain

par(mfrow=c(1,2))
plot(x$mycorrhization_percentage_HH1 ~ x$condition,ylab="% mycorrhiza",xlab="condition",main="HH1 (sensitive strain)",ylim=c(0,100))
plot(x$mycorrhization_percentage_LS4 ~ x$condition,ylab="% mycorrhiza",xlab="condition",main="LS4 (tolerant strain)",ylim=c(0,100))


lm1a <- lm(x$mycorrhization_percentage_LS4 ~ x$condition + x$mycorrhization_percentage_HH1)
lm1b <- lm(x$mycorrhization_percentage_LS4 ~ x$condition * x$mycorrhization_percentage_HH1)
lm1c <- lm(x$mycorrhization_percentage_LS4 ~ x$condition_binary + x$mycorrhization_percentage_HH1)
lm1d <- lm(x$mycorrhization_percentage_LS4 ~ x$batch)

lm2a <- lm(x$mycorrhization_percentage_HH1 ~ x$condition + x$mycorrhization_percentage_LS4)
lm2b <- lm(x$mycorrhization_percentage_HH1 ~ x$condition * x$mycorrhization_percentage_LS4)
lm2c <- lm(x$mycorrhization_percentage_HH1 ~ x$condition_binary + x$mycorrhization_percentage_LS4)
lm2d <- lm(x$mycorrhization_percentage_HH1 ~ x$batch)


# Shoot mass is not directly (significantly) dependent on (binary) condition:
lm3a  <- lm(x$shoot_mass ~ x$condition)
lm3b  <- lm(x$shoot_mass ~ x$condition_binary)

# It depends a lot on the mycorrhization rate by HH1:
# So this strain seems to provide a benefit, regardless of condition
# Or, another way of looking at it, is that if condition is difficult for plant, HH1 goes down ...
lm3c  <- lm(x$shoot_mass ~ x$mycorrhization_percentage_HH1)
lm3c  <- lm(x$shoot_mass ~ x$mycorrhization_percentage_LS4)


lm3a  <- lm(x$shoot_mass ~ x$condition_binary * x$mycorrhization_percentage_LS4 + x$mycorrhization_percentage_HH1)
lm3c  <- lm(x$shoot_mass ~ x$condition_binary * x$mycorrhization_percentage_HH1 * x$root_proportion_HH1)



# You can inspect models like so:
anova(lm1a) # etc


# Some informative plots:
plot(x$shoot_mass ~ x$condition)
plot(x$shoot_mass ~ x$condition_binary)
plot(x$root_proportion_HH1 ~ x$condition_binary)
plot(x$root_proportion_LS4 ~ x$condition_binary)

plot(x$root_logratio ~ x$condition_binary,ylab="HH1 <--  log2(LS4/HH1)   --> LS4", main = "log ratio root mass")








# Next, we will use a random forest (machine learning technique) in order to see if there is any hope predicting the condition from the data
# If not, then this means we are kind of wasting our time to try to fit a bunch of statistical models. 
accuracies  <- rep(0,1000)
importances <- rep(0,13)
for (i in 1:1000)
{
x_rf        <- x[,c("condition","shoot_mass","root_mass_HH1","root_mass_LS4","mycorrhization_percentage_HH1","mycorrhization_percentage_LS4","root_hairs_count_HH1","root_hairs_count_LS4","root_logratio","nonmycorrhized_logratio","mycorrhization_diff","root_proportion_LS4","root_proportion_HH1","shoot_proportion")]
train_index <- createDataPartition(x_rf$condition, p = runif(1,0.6,0.8), list = FALSE)
train_data  <- x_rf[train_index, ]
test_data   <- x_rf[-train_index, ]
rf_model    <- randomForest(condition ~ ., data = train_data, ntree = sample(100:300,1), importance = TRUE)

# Make predictions on the test set
predictions <- predict(rf_model, newdata = test_data)

# Evaluate performance
conf_matrix <- confusionMatrix(predictions, test_data$condition)
# print(conf_matrix)
# This gives the accuracy of the predicitons in the test data
accuracies[i] <- conf_matrix$overall[1]
importances   <- importances + importance(rf_model)[,5]
}

# Variable importance
# importance(rf_model)
# varImpPlot(rf_model)


### We will do the same, but for condition binary
accuracies_binary  <- rep(0,1000)
importances_binary <- rep(0,13)
for (i in 1:1000)
{  
x_rf        <- x[,c("condition_binary","shoot_mass","root_mass_HH1","root_mass_LS4","mycorrhization_percentage_HH1","mycorrhization_percentage_LS4","root_hairs_count_HH1","root_hairs_count_LS4","root_logratio","nonmycorrhized_logratio","mycorrhization_diff","root_proportion_LS4","root_proportion_HH1","shoot_proportion")]
train_index <- createDataPartition(x_rf$condition_binary, p = runif(1,0.6,0.8), list = FALSE)
train_data  <- x_rf[train_index, ]
test_data   <- x_rf[-train_index, ]

rf_model    <- randomForest(condition_binary ~ ., data = train_data, ntree = sample(100:300,1), importance = TRUE)

# Make predictions on the test set
predictions <- predict(rf_model, newdata = test_data)

# Evaluate performance
conf_matrix <- confusionMatrix(predictions, test_data$condition)
#print(conf_matrix)
accuracies_binary[i] <- conf_matrix$overall[1]
importances_binary   <- importances_binary + importance(rf_model)[,4]

}

# Density of accuracies. Cleary not accurate.
plot(density(accuracies,bw=0.03),xlim=c(0,1),ylim=c(0,4), xlab="accuracy of classification",main="accuracies of classification by random forest")
lines(density(accuracies_binary,bw=0.03),col="red")

names(importances)        <- names(importance(rf_model)[,4])
names(importances_binary) <- names(importance(rf_model)[,4])


# Variable importance
# importance(rf_model)
# varImpPlot(rf_model)


# CONCLUSION FOR RANDOM FORESTS:
# Classification into conditions is not possible
# Classification into binary conditions (zinc vs no zinc) is better (around 75%), but also not fantastic. 
# Some conclusions from importances (lower values are higher ranks):
# * The most important variables are related to shoot mass and shoot/root proportion
# * The root hair counts add very little information
# * Mycorrhization rates are also not super informative, especially not for the HH1 side











# SEM
mod_shoot <- lm(shoot_mass ~ 
                  condition_binary * mycorrhization_percentage_HH1 * mycorrhization_percentage_LS4, 
                  data = x_rf)

mod_HH1 <- lm(
  mycorrhization_percentage_HH1 ~ condition_binary + mycorrhization_percentage_LS4,
  data = x_rf
)

mod_LS4 <- lm(
  mycorrhization_percentage_LS4 ~ condition_binary + mycorrhization_percentage_HH1,
  data = x_rf
)

full_sem <- psem(
  mod_shoot,
  mod_HH1,
  mod_LS4
)

summary(full_sem)

