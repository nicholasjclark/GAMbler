# Simple exploration of MyHospitals ED data
library(readaihw) # pak::pkg_install("RWParsons/readaihw")
library(tidyverse)
library(viridis)
library(scales)
library(leaflet)
library(sf)

# Get the data
hosp_locs <- get_hospital_mappings()
ed_presentations <- read_flat_data_extract("MYH-ED")

# Clean and prepare the data - Focus on NSW and ACT
ed_clean <- ed_presentations %>%
  mutate(
    value = as.numeric(value),
    date = as.Date(reporting_end_date),
    triage_category = reported_measure_category_name,
    state = mapped_state
  ) %>%
  # Include NSW and ACT
  filter(!is.na(value), value > 0, state %in% c("NSW", "ACT"))

# Plot 1: NSW + ACT ED Presentations Over Time
nsw_ed_trend <- ed_clean %>% 
  filter(reporting_unit_type_name == "State") %>% 
  group_by(date) %>% 
  summarise(total_presentations = sum(value, na.rm = TRUE), .groups = 'drop')

ggplot(nsw_ed_trend, aes(x = date, y = total_presentations)) +
  geom_line(size = 1.2, color = "darkred") +
  geom_point(size = 2, color = "darkred") +
  labs(title = "NSW + ACT ED Presentations Over Time",
       x = "Date", 
       y = "Total ED Presentations") +
  theme_classic() +
  scale_y_continuous(labels = comma_format())

# Plot 2: NSW + ACT Triage Category Distribution
triage_summary <- ed_clean %>%
  group_by(triage_category) %>%
  summarise(total_presentations = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(total_presentations))

ggplot(triage_summary, aes(x = reorder(triage_category, total_presentations), 
                                 y = total_presentations, fill = triage_category)) +
  geom_col() +
  coord_flip() +
  labs(title = "NSW + ACT ED Presentations by Triage Category",
       x = "Triage Category", 
       y = "Total Presentations") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = comma_format())


# Plot 3: Interactive Map of NSW + ACT Hospitals with ED Activity
nsw_hospital_activity <- ed_clean %>%
  filter(reporting_unit_type_name == "Hospital") %>%
  group_by(reporting_unit_code, reporting_unit_name) %>%
  summarise(
    total_presentations = sum(value, na.rm = TRUE),
    avg_annual = mean(value, na.rm = TRUE),
    years_reporting = n_distinct(year(date)),
    .groups = 'drop'
  ) %>%
  left_join(
    hosp_locs %>% 
      filter(state %in% c("New South Wales", "Australian Capital Territory")) %>%
      mutate(
        lat = as.numeric(latitude),
        lon = as.numeric(longitude)
      ) %>%
      select(code, name, lat, lon, sector, local_hospital_network_lhn),
    by = c("reporting_unit_code" = "code")
  ) %>%
  filter(!is.na(lat), !is.na(lon), total_presentations > 0)

# Load Australian state boundaries for mask layer
aus_states <- read_sf("https://raw.githubusercontent.com/rowanhogan/australian-states/master/states.geojson")

# Create interactive map with minimal basemap
map <- leaflet(nsw_hospital_activity) %>%
  # Add state boundaries with NSW and ACT highlighted
  addPolygons(
    data = aus_states,
    fillColor = ~ifelse(STATE_NAME %in% c("New South Wales", "Australian Capital Territory"), "transparent", "black"),
    fillOpacity = ~ifelse(STATE_NAME %in% c("New South Wales", "Australian Capital Territory"), 0, 0.5),
    color = "white",
    weight = 0.5,
    opacity = 0.5
  ) %>%
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    # Scale by annual average
    radius = ~sqrt(total_presentations / years_reporting) / 18,
    color = "white",
    fillColor = "darkred",
    fillOpacity = 0.4,
    stroke = TRUE,
    weight = 1,
    popup = ~paste(
      "<b>", name, "</b><br/>",
      "Total ED Presentations: ", scales::comma(total_presentations), "<br/>",
      "Network: ", local_hospital_network_lhn, "<br/>",
      "Years Reporting: ", years_reporting
    )
  ) %>%
  addProviderTiles("CartoDB.Positron") %>%  # Clean, minimal basemap
  # Add combined title and legend
  addControl(
    html = '<div style="background: white; padding: 8px; border-radius: 5px; box-shadow: 0 0 10px rgba(0,0,0,0.2); font-family: Arial, sans-serif;">
             <div style="font-size: 12px; font-weight: bold; color: #333; margin-bottom: 10px; border-bottom: 1px solid #ddd; padding-bottom: 8px;">
               NSW & ACT Emergencies
             </div>
             <div style="margin: 5px 0; display: flex; align-items: center;">
               <div style="width: 20px; height: 10px; display: flex; justify-content: center; align-items: center; margin-right: 8px;">
                 <span style="width: 20px; height: 20px; border-radius: 50%; background: darkred; display: block;"></span>
               </div>
               <span style="font-size: 11px;">High Volume</span>
             </div>
             <div style="margin: 8px 0; display: flex; align-items: center;">
               <div style="width: 20px; height: 10px; display: flex; justify-content: center; align-items: center; margin-right: 8px;">
                 <span style="width: 12px; height: 12px; border-radius: 50%; background: darkred; display: block;"></span>
               </div>
               <span style="font-size: 11px;">Medium Volume</span>
             </div>
             <div style="margin: 8px 0; display: flex; align-items: center;">
               <div style="width: 20px; height: 10px; display: flex; justify-content: center; align-items: center; margin-right: 8px;">
                 <span style="width: 6px; height: 6px; border-radius: 50%; background: darkred; display: block;"></span>
               </div>
               <span style="font-size: 11px;">Low Volume</span>
             </div>
           </div>',
    position = "topright"
  ) %>%
  # Center on Dubbo with fixed zoom
  setView(lng = 148.6019, lat = -32.2431, zoom = 6)

print(map)


# TO DO
# 1. Add Hospital Network polygons as a toggle layer (first research where to find these)
# 2. Any way to have place names show ABOVE points for readability?
# 3. Add toggle switch or dropdown to modify map: 
#     - compute trends over time, relative to population (get population layer from ABS, look her: https://digital.atlas.gov.au/maps/digitalatlas::abs-australian-population-grid-2024/about) and colour points by slope of trend
#     - replace points with spatial smooth predictions from mgcv using soap film smoother
# 4. Add histogram of total emergencies in last year
#     - when user clicks a point, a vertical line is added to histogram at its value