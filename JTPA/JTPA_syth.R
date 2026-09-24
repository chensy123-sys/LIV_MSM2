rm(list = ls())
source('LIV_xgboost.R'); source("Script/par/par_overall_continuous_jtpa.R")
# source('LIV_mgcv.R')
set.seed(2025)
my_data_long <- sim_LIV_overall(n = 10000,D)
# fit_LIV_00_local <- DML_LIV(mydat = my_data_long, action_list = c(0, 0), estimand = "local")
# fit_LIV_01_local <- DML_LIV(mydat = my_data_long, action_list = c(0, 1), estimand = "local")
# fit_LIV_10_local <- DML_LIV(mydat = my_data_long, action_list = c(1, 0), estimand = "local")
# fit_LIV_11_local <- DML_LIV(mydat = my_data_long, action_list = c(1, 1), estimand = "local")

fit_LIV_00_overall <- DML_LIV(mydat = my_data_long, action_list = c(0, 0), estimand = "overall")
fit_LIV_01_overall <- DML_LIV(mydat = my_data_long, action_list = c(0, 1), estimand = "overall")
fit_LIV_10_overall <- DML_LIV(mydat = my_data_long, action_list = c(1, 0), estimand = "overall")
fit_LIV_11_overall <- DML_LIV(mydat = my_data_long, action_list = c(1, 1), estimand = "overall")


# fit_LIV_1p_local <- DML_LIV(mydat = my_data_long, action_list = list(
#   function(fY) return(1),
#   function(fY) {return(as.numeric(fY$L<0.8368623))}), estimand = "local")
# fit_LIV_1n_local <- DML_LIV(mydat = my_data_long, action_list = list(
#   function(fY) return(1),
#   function(fY) {return(as.numeric(fY$L>=0.8368623))}), estimand = "local")
# fit_LIV_0p_local <- DML_LIV(mydat = my_data_long, action_list = list(
#   function(fY) return(0),
#   function(fY) {return(as.numeric(fY$L<0.8368623))}), estimand = "local")
# fit_LIV_0n_local <- DML_LIV(mydat = my_data_long, action_list = list(
#   function(fY) return(0),
#   function(fY) {return(as.numeric(fY$L<=0.8368623))}), estimand = "local")
# 
# 
fit_LIV_1p_overall <- DML_LIV(mydat = my_data_long, action_list = list(
  function(fY) return(1),
  function(fY) {return(as.numeric(fY$L<0.8368623))}), estimand = "overall")
fit_LIV_1n_overall <- DML_LIV(mydat = my_data_long, action_list = list(
  function(fY) return(1),
  function(fY) {return(as.numeric(fY$L>=0.8368623))}), estimand = "overall")
fit_LIV_0p_overall <- DML_LIV(mydat = my_data_long, action_list = list(
  function(fY) return(0),
  function(fY) {return(as.numeric(fY$L<0.8368623))}), estimand = "overall")
fit_LIV_0n_overall <- DML_LIV(mydat = my_data_long, action_list = list(
  function(fY) return(0),
  function(fY) {return(as.numeric(fY$L>=0.8368623))}), estimand = "overall")


# pred_LIV_local <- rbind(
  # est_00 = predict_LIV(fit_LIV_00_local) %>% select(-lower,-upper),
  # est_01 = predict_LIV(fit_LIV_01_local) %>% select(-lower,-upper),
  # est_0p = predict_LIV(fit_LIV_0p_local) %>% select(-lower,-upper),
  # est_0n = predict_LIV(fit_LIV_0n_local) %>% select(-lower,-upper),
  # est_10 = predict_LIV(fit_LIV_10_local) %>% select(-lower,-upper),
  # est_11 = predict_LIV(fit_LIV_11_local) %>% select(-lower,-upper)
  # est_1p = predict_LIV(fit_LIV_1p_local) %>% select(-lower,-upper),
  # est_1n = predict_LIV(fit_LIV_1n_local) %>% select(-lower,-upper)
# )
pred_LIV_overall <- rbind(
  est_00 = predict_LIV(fit_LIV_00_overall) %>% select(-lower,-upper),
  est_01 = predict_LIV(fit_LIV_01_overall) %>% select(-lower,-upper),
  # est_0p = predict_LIV(fit_LIV_0p_overall) %>% select(-lower,-upper),
  # est_0n = predict_LIV(fit_LIV_0n_overall) %>% select(-lower,-upper),
  est_10 = predict_LIV(fit_LIV_10_overall) %>% select(-lower,-upper),
  est_11 = predict_LIV(fit_LIV_11_overall) %>% select(-lower,-upper)
  # est_1p = predict_LIV(fit_LIV_1p_overall) %>% select(-lower,-upper),
  # est_1n = predict_LIV(fit_LIV_1n_overall) %>% select(-lower,-upper)
)
rbind(
  pred_LIV_overall %>% filter(time==0) %>% select (-time),
  (pred_LIV_overall %>% filter(time==1))[1:4,]%>% select (-time) 
) %>% t() %>% data.frame() %>%
  mutate(across(everything(), ~ round(.x, 4)))


# rbind(
#   pred_LIV_local %>% filter(time==0) %>% select (-time),
#   (pred_LIV_local %>% filter(time==1))[1:4,]%>% select (-time) 
# ) %>% t() %>% data.frame() %>%
#   mutate(across(everything(), ~ round(.x, 4)))
