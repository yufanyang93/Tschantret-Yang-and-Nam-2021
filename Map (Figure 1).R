library(ggmap)
library(ggplot2)
library(haven)
library(dplyr)
library(tidyr)
library(raster)
library(tmap)
library(sf)
library(spData)
library(rgdal)
library(RColorBrewer)
library(countrycode)
library(tigris)
library(rnaturalearth)
library(rnaturalearthdata)

rm(list = ls())

##Pre 1945 Map##
terr_group<-read_dta("Historical_Terrorism_Data.dta")

terr_group$iso_a2<-countrycode(terr_group$country, "country.name", "iso2c")
terr_group$iso_a2<-ifelse(terr_group$country=="Syria/Turkey", "TR", terr_group$iso_a2)
terr_group$iso_a2<-ifelse(terr_group$country=="Algiers", "DZ", terr_group$iso_a2)
terr_group$iso_a2<-ifelse(terr_group$country=="Columbia", "CO", terr_group$iso_a2)
terr_group$iso_a2<-ifelse(terr_group$country=="Uraguay", "UY", terr_group$iso_a2)
terr_group$iso_a2<-ifelse(terr_group$country=="Wales", "GB", terr_group$iso_a2)

terr_group<-terr_group%>%group_by(iso_a2)%>%mutate(countrytotal=sum(terronset))
terr_group<-terr_group%>%group_by(city)%>%mutate(citytotal=sum(terronset))
terr_group_point<-terr_group%>%filter(terronset==1)
terr_group_ymca<-terr_group%>%group_by(city)%>%
  mutate(ymca2=ifelse(sum(ymca>=1), 1, 0))%>%
  distinct(city, longitude, latitude, ymca2, .keep_all = TRUE)%>%
  filter(ymca2==1)

wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
wmap_df<-fortify(wmap)
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="#e6e8ed"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))

countries <- readOGR("ne_110m_admin_0_countries", layer="ne_110m_admin_0_countries") 
countries_df<-fortify(countries)

grat <- readOGR("ne_110m_graticules_all", layer="ne_110m_graticules_15") 
grat_df <- fortify(grat)

bbox <- readOGR("ne_110m_graticules_all", layer="ne_110m_wgs84_bounding_box") 
bbox_df<- fortify(bbox)

polity<-read_dta("p5.dta")

polity$polity2<-ifelse(polity$polity2<=-11, NA, polity$polity2)
polity$ISO_A2<-countrycode(polity$country, "country.name", "iso2c")

polity$ISO_A2<-ifelse(polity$country=="Bavaria", "DE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Czechoslovakia", "CZ", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Germany East", "DE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="South Vietnam", "VN", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yugoslavia", "YU", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yemen North", "YE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yemen South", "YE", polity$ISO_A2)

polity<-polity%>%filter(year>=1860 & year<=1945)%>%
  group_by(ISO_A2)%>%filter(!is.na(polity2))%>%mutate(ave_polity=mean(polity2))

polity<-polity%>%distinct(country, ISO_A2, ave_polity)

world_polity_df<-fortify(countries, region = "ISO_A2")
world_polity_df<-left_join(world_polity_df, polity, by=c('id'='ISO_A2'))
world_polity_df<-world_polity_df%>%group_by(id)%>%mutate(regime=ifelse(ave_polity>=6, "Democracy", ifelse(ave_polity %in% -5:5, "Anocracy", "Autocracy")))

worldmap<-ggplot(wmap_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=countries_df, aes(long,lat, group=group)) + 
  geom_polygon(data = world_polity_df, aes(long, lat, group=group, fill=regime))+
  geom_path(data=countries_df, aes(long,lat, group=group), color="white", size=0.3) +
  geom_path(data=grat_df, aes(long, lat, group=group, fill=NULL), size=0.1, linetype="dashed", color="grey50") +
  geom_point(data=terr_group_point, aes(longitude, latitude, group=NULL, fill=NULL, size=citytotal), color="#3655B3", alpha=I(7/10), show.legend = FALSE) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("#EA9A62", "#A45445", "#6BB300"), na.value = "#A6A1A2")+
  scale_size_continuous(range=c(1, 8), guide="none")+
  labs(title="Historical Terrorist Groups and Regime Types Pre-1945", fill="Regime")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("pre1945.png", width=12.5, height=8.25, dpi=300)

##Post 1945 Map##
rm(list = ls())

polity<-read_dta("p5.dta")

polity$polity2<-ifelse(polity$polity2<=-11, NA, polity$polity2)
polity$ISO_A2<-countrycode(polity$country, "country.name", "iso2c")

polity$ISO_A2<-ifelse(polity$country=="Bavaria", "DE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Czechoslovakia", "CZ", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Germany East", "DE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="South Vietnam", "VN", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yugoslavia", "YU", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yemen North", "YE", polity$ISO_A2)
polity$ISO_A2<-ifelse(polity$country=="Yemen South", "YE", polity$ISO_A2)

polity<-polity%>%filter(year>=1946 & year<=1970)%>%
  group_by(ISO_A2)%>%filter(!is.na(polity2))%>%mutate(ave_polity=mean(polity2))

polity<-polity%>%distinct(country, ISO_A2, ave_polity)

wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
wmap_df<-fortify(wmap)
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="#e6e8ed"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))

countries <- readOGR("ne_110m_admin_0_countries", layer="ne_110m_admin_0_countries") 
countries_df<-fortify(countries)

grat <- readOGR("ne_110m_graticules_all", layer="ne_110m_graticules_15") 
grat_df <- fortify(grat)

bbox <- readOGR("ne_110m_graticules_all", layer="ne_110m_wgs84_bounding_box") 
bbox_df<- fortify(bbox)

world_polity_df<-fortify(countries, region = "ISO_A2")
world_polity_df<-left_join(world_polity_df, polity, by=c('id'='ISO_A2'))
world_polity_df<-world_polity_df%>%group_by(id)%>%mutate(regime=ifelse(ave_polity>=6, "Democracy", ifelse(ave_polity %in% -5:5, "Anocracy", "Autocracy")))

terr_group<-read.csv("terrorist-groups.csv")
terr_group<-terr_group%>%filter(StrYear>=1946)%>%unite(location, Country, City, sep = ", ", remove = FALSE)

register_google(key = "MYKEY")
terr_group<-terr_group%>%mutate_geocode(location, output = c("more"), source = c("google", "dsk"))
terr_group_point<-terr_group%>%filter(!is.na(City))%>%group_by(City)%>%mutate(citytotal=n())

worldmap<-ggplot(wmap_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=countries_df, aes(long,lat, group=group)) + 
  geom_polygon(data = world_polity_df, aes(long, lat, group=group, fill=regime))+
  geom_path(data=countries_df, aes(long,lat, group=group), color="white", size=0.3) +
  geom_path(data=grat_df, aes(long, lat, group=group, fill=NULL), size=0.1, linetype="dashed", color="grey50") +
  geom_point(data=terr_group_point, aes(lon, lat, group=NULL, fill=NULL, size=citytotal), color="#3655B3", alpha=I(8/10), show.legend = FALSE) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("#EA9A62", "#A45445", "#6BB300"), na.value = "#A6A1A2")+
  scale_size_continuous(range=c(1, 8), guide="none")+
  labs(title="Historical Terrorist Groups and Regime Types Post-1945", fill="Regime")+
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("post1945.png", width=12.5, height=8.25, dpi=300)
