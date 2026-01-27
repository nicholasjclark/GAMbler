# From Academia to Industry: Building Production Data Pipelines with R, SQL, 
# and Modern Cloud Tools
# Many data scientists think they need to abandon R for Python to work with 
# modern data warehouses, but here I demonstrate how to build a complete 
# production pipeline—extracting 50,000 NYC taxi trips directly from cloud 
# storage, engineering features with SQL window functions, training 
# predictive models in R, and deploying results back to a data warehouse—all 
# in under 100 lines of code. This pipeline showcases the exact skills tech 
# companies seek: efficient data extraction without downloading massive files, 
# complex SQL transformations using CTEs and analytical functions, statistical 
# modeling integrated with database operations, and bidirectional data flow 
# between analytical and production systems. Drawing from my experience 
# developing the mvgam package for ecological forecasting (processing millions 
#                                                          of observations 
#                                                          across hierarchical 
#                                                          spatiotemporal 
#                                                          structures), 
# I show how the same rigorous statistical approaches that power academic 
# research translate directly to industry-scale data engineering. Whether 
# you're analyzing ecological time series or predicting taxi fares, the 
# core challenge remains the same: efficiently transforming raw data into 
# actionable insights at scale—and R, combined with DuckDB's SQL engine, 
# proves to be an unexpectedly powerful solution for both.

# NYC Taxi Pipeline - Efficient Version (No Full Download)
library(duckdb)
library(DBI)
library(tidyverse)

# 1. SETUP: Connection with httpfs
con <- dbConnect(duckdb())
dbExecute(con, "INSTALL httpfs; LOAD httpfs")

# 2. Query remote file
sample_data <- dbGetQuery(con, "
  SELECT * FROM read_parquet('https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2024-01.parquet')
  WHERE trip_distance > 0 
    AND fare_amount > 0
    AND DATE(tpep_pickup_datetime) = '2024-01-15'  -- Just one day!
  LIMIT 50000
")

# Write sample to local table for faster processing
dbWriteTable(con, "taxi_sample", sample_data, overwrite = TRUE)

# 3. EXPLORATION: Quick overview
overview <- dbGetQuery(con, "
  SELECT 
    COUNT(*) as trips,
    DATE(MIN(tpep_pickup_datetime)) as first_date,
    DATE(MAX(tpep_pickup_datetime)) as last_date,
    ROUND(AVG(trip_distance), 2) as avg_distance,
    ROUND(AVG(total_amount), 2) as avg_fare
  FROM taxi_sample
")
print(overview)

# 4. FEATURE ENGINEERING with SQL
dbExecute(con, "
  CREATE OR REPLACE TABLE taxi_features AS
  SELECT 
    tpep_pickup_datetime as pickup_time,
    trip_distance,
    fare_amount,
    tip_amount,
    PULocationID as pickup_zone,
    -- Calculate derived features
    EXTRACT(HOUR FROM tpep_pickup_datetime) as hour,
    EXTRACT(DOW FROM tpep_pickup_datetime) as day_of_week,
    ROUND(fare_amount / NULLIF(trip_distance, 0), 2) as fare_per_mile,
    ROUND(tip_amount / NULLIF(fare_amount, 0) * 100, 1) as tip_pct,
    CASE 
      WHEN EXTRACT(HOUR FROM tpep_pickup_datetime) BETWEEN 7 AND 10 THEN 'Morning'
      WHEN EXTRACT(HOUR FROM tpep_pickup_datetime) BETWEEN 17 AND 20 THEN 'Evening'
      ELSE 'Off-Peak'
    END as rush_hour
  FROM taxi_sample
  WHERE trip_distance BETWEEN 0.5 AND 50
    AND fare_amount BETWEEN 5 AND 200
")

# 5. SQL ANALYTICS: Zone performance
zone_analysis <- dbGetQuery(con, "
  WITH zone_metrics AS (
    SELECT 
      pickup_zone,
      rush_hour,
      COUNT(*) as trips,
      ROUND(AVG(fare_per_mile), 2) as avg_fare_per_mile,
      ROUND(AVG(tip_pct), 1) as avg_tip_pct
    FROM taxi_features
    GROUP BY pickup_zone, rush_hour
    HAVING COUNT(*) >= 10
  )
  SELECT 
    pickup_zone,
    rush_hour,
    trips,
    avg_fare_per_mile,
    avg_tip_pct,
    RANK() OVER (PARTITION BY rush_hour ORDER BY trips DESC) as zone_rank
  FROM zone_metrics
  ORDER BY rush_hour, zone_rank
  LIMIT 15
")

print(zone_analysis)

# 6. R MODELING: Quick tip prediction model
model_data <- dbGetQuery(con, "
  SELECT * FROM taxi_features 
  WHERE tip_pct IS NOT NULL
")

# Simple model
tip_model <- lm(tip_pct ~ hour + day_of_week + trip_distance + fare_per_mile + rush_hour,
                data = model_data)

# 7. PREDICTIONS
model_data$predicted_tip <- predict(tip_model)
model_data$error <- abs(model_data$tip_pct - model_data$predicted_tip)

# 8. WRITE BACK
predictions <- model_data %>%
  select(pickup_time, tip_pct, predicted_tip, error) %>%
  slice_head(n = 1000)

dbWriteTable(con, "predictions", predictions, overwrite = TRUE)

# 9. PERFORMANCE SUMMARY
performance <- dbGetQuery(con, "
  SELECT 
    COUNT(*) as n,
    ROUND(AVG(error), 2) as mae,
    ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY error), 2) as median_error
  FROM predictions
")

cat("\n=== RESULTS ===\n")
cat("Samples processed:", nrow(model_data), "\n")
cat("Model R²:", round(summary(tip_model)$r.squared, 3), "\n")
cat("Mean Absolute Error:", performance$mae, "%\n")

dbDisconnect(con)
