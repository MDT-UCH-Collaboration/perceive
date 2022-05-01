library(tidymodels)
library(modeltime)
library(tidyverse)
library(timetk)

# Used to convert plots from interactive to static
interactive = TRUE

# Read data
bike_transactions_tbl <- bike_sharing_daily %>%
  select(dteday, cnt) %>%
  set_names(c("date", "value")) 

bike_transactions_tbl

bike_transactions_tbl %>%
  plot_time_series(date, value, .interactive = interactive)


splits <- bike_transactions_tbl %>%
  time_series_split(assess = "3 months", cumulative = TRUE)


splits %>%
  tk_time_series_cv_plan() %>%
  plot_time_series_cv_plan(date, value, .interactive = interactive)
