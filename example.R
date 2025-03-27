# load packages
source("cdli_aces_code.R")

# save seed for reproducibility
set.seed(1)

# ================================== #
#           Simulate Data            # 
# ================================== #

data <- aces:::cobb_douglas_XnY1 (
  N = 300,
  nX = 3
)

# simulate three different countries
data[, "country"] <- sample (
  c("country A", "country B", "country C"), 
  size = 300, 
  replace = TRUE
  )

# simulate three different years
data[, "year"] <- sample (
  c("year 1", "year 2", "year 3"), 
  size = 300, 
  replace = TRUE
)

# ================================== #
#        Build Reference Group       # 
# ================================== #

# margin error
e <- 5 / 100

# confidence level
z <- qnorm(1 - (0.05 / 2))

# proportion
p <- 0.5

# population size
N <- nrow(data)

# sample size
n <- ceiling((N * z ^ 2 * p * (1 - p)) / (e ^ 2 * (N - 1) + z ^ 2 * p * (1 - p)))

# proportions in data
refer_table <- ceiling(prop.table(table(data$country, data$year)) * n)

# number of observations for each country and year
count_table <- table(data$country, data$year)

# initialize reference data.frame
group_R <- data %>%
  filter(country == "country A", year == "year 1") %>%
  sample_n(15)

# list of countries
countries <- c("country A", "country B", "country C")

# list of years
years <- c("year 1", "year 2", "year 3")

for (i in 1:3) {
  for (j in 1:3) {
    
    if (i == 1 && j == 1) next
    
    # sample size for "country-year" combination
    sample_size <- refer_table[i, j]
    
    # sample with the required length
    sample_data <- data %>%
      filter(country == countries[i], year == years[j]) %>%
      sample_n(sample_size)
    
    # add sample_data to group_R
    group_R <- rbind(group_R, sample_data)
    
  }
}

# ================================== #
#   Reference Group: Set Technology  #
# ================================== #

# build and ACES model with the reference group
aces_model_R <- aces (
  data = group_R,
  x = 1:3,
  y = 4,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "inter_cost" = 1
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  max_terms = nrow(group_R),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 1
)

group_R_aces <- cbind (
  aces_model_R[["technology"]][["aces"]][["xmat"]],
  aces_model_R[["technology"]][["aces"]][["ymat"]]
)

# ================================== #
# Camanho-Dyson Luenberger Indicator #
# ================================== #

# We show an example of how to calculate the Camanho-Dyson Luenberger Indicator
# for Country A and Country B between Year 1 and Year 2

# ==
# Filter data: year 1
# ==

# Country A for year 1
group_A_t1 <- data %>%
  filter(country == "country A", year == "year 1") %>%
  sample_n(min(count_table[, "year 1"]))

# Fit the model
aces_model_A_t1 <- aces (
  data = group_A_t1,
  x = 1:3,
  y = 4,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "inter_cost" = 1
    ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
    ),
  max_terms = nrow(group_A_t1),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 1
  )

# Set technology
group_A_t1_aces <- cbind (
  aces_model_A_t1[["technology"]][["aces"]][["xmat"]],
  aces_model_A_t1[["technology"]][["aces"]][["ymat"]]
  )

# Country B for year 1
group_B_t1 <- data %>%
  filter(country == "country B", year == "year 1") %>%
  sample_n(min(count_table[, "year 1"]))

# Fit the model
aces_model_B_t1 <- aces (
  data = group_B_t1,
  x = 1:3,
  y = 4,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "inter_cost" = 1
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  max_terms = nrow(group_B_t1),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 1
  )

# Set technology
group_B_t1_aces <- cbind (
  aces_model_B_t1[["technology"]][["aces"]][["xmat"]],
  aces_model_B_t1[["technology"]][["aces"]][["ymat"]]
  )

# ==
# Camano Dyson Luenberger Indicator: year 1
# ==

ACES_cdli_t1 <- cdli (
  tech_xmat_R = as.matrix(group_R_aces[, 1:3]),
  tech_ymat_R = as.matrix(group_R_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t1[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t1[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t1[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t1[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

# ==
# Filter data: year 2
# ==

# Country A for year 2
group_A_t2 <- data %>%
  filter(country == "country A", year == "year 2") %>%
  sample_n(min(count_table[, "year 2"]))

# Fit the model
aces_model_A_t2 <- aces (
  data = group_A_t2,
  x = 1:3,
  y = 5,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "inter_cost" = 1
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  max_terms = nrow(group_A_t2),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 1
)

# Set technology
group_A_t2_aces <- cbind (
  aces_model_A_t2[["technology"]][["aces"]][["xmat"]],
  aces_model_A_t2[["technology"]][["aces"]][["ymat"]]
)

# Country B for year 2
group_B_t2 <- data %>%
  filter(country == "country B", year == "year 2") %>%
  sample_n(min(count_table[, "year 2"]))

# Fit the model
aces_model_B_t2 <- aces (
  data = group_B_t2,
  x = 1:3,
  y = 5,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "inter_cost" = 1
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  max_terms = nrow(group_B_t2),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 1
  )

# Set technology
group_B_t2_aces <- cbind (
  aces_model_B_t2[["technology"]][["aces"]][["xmat"]],
  aces_model_B_t2[["technology"]][["aces"]][["ymat"]]
)

# ==
# Camano Dyson Luenberger Indicator: year 2
# ==

ACES_cdli_t2 <- cdli (
  tech_xmat_R = as.matrix(group_R_aces[, 1:3]),
  tech_ymat_R = as.matrix(group_R_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t2[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t2[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t2[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t2[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

# ==
# Pseudo-Panel Luenberger Indicator
# ==

ppli_ACES <- ACES_cdli_t1 - ACES_cdli_t2

# ==
# Efficiency Gap
# ==

egap_ppli_ACES <- egap_ppli (
  tech_xmat_A_t1 = as.matrix(group_A_t1_aces[, 1:3]),
  tech_ymat_A_t1 = as.matrix(group_A_t1_aces[, 4]),
  tech_xmat_B_t1 = as.matrix(group_B_t1_aces[, 1:3]),
  tech_ymat_B_t1 = as.matrix(group_B_t1_aces[, 4]),
  tech_xmat_A_t2 = as.matrix(group_A_t2_aces[, 1:3]),
  tech_ymat_A_t2 = as.matrix(group_A_t2_aces[, 4]),
  tech_xmat_B_t2 = as.matrix(group_B_t2_aces[, 1:3]),
  tech_ymat_B_t2 = as.matrix(group_B_t2_aces[, 4]),
  eval_xmat_A_t1 = as.matrix(group_A_t1[, 1:3]),
  eval_ymat_A_t1 = as.matrix(group_A_t1[, c(4)]),
  eval_xmat_B_t1 = as.matrix(group_B_t1[, 1:3]),
  eval_ymat_B_t1 = as.matrix(group_B_t1[, c(4)]),
  eval_xmat_A_t2 = as.matrix(group_A_t2[, 1:3]),
  eval_ymat_A_t2 = as.matrix(group_A_t2[, c(4)]),
  eval_xmat_B_t2 = as.matrix(group_B_t2[, 1:3]),
  eval_ymat_B_t2 = as.matrix(group_B_t2[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
  )

egap_cdli_t1 <- egap_cdli (
  tech_xmat_A_t = as.matrix(group_A_t1_aces[, 1:3]),
  tech_ymat_A_t = as.matrix(group_A_t1_aces[, 4]),
  tech_xmat_B_t = as.matrix(group_B_t1_aces[, 1:3]),
  tech_ymat_B_t = as.matrix(group_B_t1_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t1[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t1[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t1[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t1[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

egap_cdli_t1_ACES <- egap_cdli_t1[["EG"]]

egap_cdli_t2 <- egap_cdli (
  tech_xmat_A_t = as.matrix(group_A_t2_aces[, 1:3]),
  tech_ymat_A_t = as.matrix(group_A_t2_aces[, 4]),
  tech_xmat_B_t = as.matrix(group_B_t2_aces[, 1:3]),
  tech_ymat_B_t = as.matrix(group_B_t2_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t2[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t2[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t2[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t2[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

egap_cdli_t2_ACES <- egap_cdli_t2[["EG"]]

# ==
# Technological gap
# ==

tgap_ppli_ACES <- tgap_ppli (
  tech_xmat_R = as.matrix(group_R_aces[, 1:3]),
  tech_ymat_R = as.matrix(group_R_aces[, 4]),
  tech_xmat_A_t1 = as.matrix(group_A_t1_aces[, 1:3]),
  tech_ymat_A_t1 = as.matrix(group_A_t1_aces[, 4]),
  tech_xmat_B_t1 = as.matrix(group_B_t1_aces[, 1:3]),
  tech_ymat_B_t1 = as.matrix(group_B_t1_aces[, 4]),
  tech_xmat_A_t2 = as.matrix(group_A_t2_aces[, 1:3]),
  tech_ymat_A_t2 = as.matrix(group_A_t2_aces[, 4]),
  tech_xmat_B_t2 = as.matrix(group_B_t2_aces[, 1:3]),
  tech_ymat_B_t2 = as.matrix(group_B_t2_aces[, 4]),
  eval_xmat_A_t1 = as.matrix(group_A_t1[, 1:3]),
  eval_ymat_A_t1 = as.matrix(group_A_t1[, c(4)]),
  eval_xmat_B_t1 = as.matrix(group_B_t1[, 1:3]),
  eval_ymat_B_t1 = as.matrix(group_B_t1[, c(4)]),
  eval_xmat_A_t2 = as.matrix(group_A_t2[, 1:3]),
  eval_ymat_A_t2 = as.matrix(group_A_t2[, c(4)]),
  eval_xmat_B_t2 = as.matrix(group_B_t2[, 1:3]),
  eval_ymat_B_t2 = as.matrix(group_B_t2[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

tgap_cdli_t1 <- tgap_cdli (
  tech_xmat_R = as.matrix(group_R_aces[, 1:3]),
  tech_ymat_R = as.matrix(group_R_aces[, 4]),
  tech_xmat_A_t = as.matrix(group_A_t1_aces[, 1:3]),
  tech_ymat_A_t = as.matrix(group_A_t1_aces[, 4]),
  tech_xmat_B_t = as.matrix(group_B_t1_aces[, 1:3]),
  tech_ymat_B_t = as.matrix(group_B_t1_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t1[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t1[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t1[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t1[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

tgap_cdli_t1_ACES <- tgap_cdli_t1[["tg"]]

tgap_cdli_t2 <- tgap_cdli (
  tech_xmat_R = as.matrix(group_R_aces[, 1:3]),
  tech_ymat_R = as.matrix(group_R_aces[, 4]),
  tech_xmat_A_t = as.matrix(group_A_t2_aces[, 1:3]),
  tech_ymat_A_t = as.matrix(group_A_t2_aces[, 4]),
  tech_xmat_B_t = as.matrix(group_B_t2_aces[, 1:3]),
  tech_ymat_B_t = as.matrix(group_B_t2_aces[, 4]),
  eval_xmat_A_t = as.matrix(group_A_t2[, 1:3]),
  eval_ymat_A_t = as.matrix(group_A_t2[, c(4)]),
  eval_xmat_B_t = as.matrix(group_B_t2[, 1:3]),
  eval_ymat_B_t = as.matrix(group_B_t2[, c(4)]),
  Gx = colMeans(data[, 1:3]),
  Gy = mean(data[, 4]),
  convexity = TRUE,
  returns = "variable"
)

tgap_cdli_t2_ACES <- tgap_cdli_t2[["tg"]]