library(tidyverse)
library(sf)

leido <- read.csv("base_madre_parcelaLRP.csv", dec = ",") 

datos <- leido %>% mutate(fecha = as.POSIXct(fecha,format = "%d/%m/%Y")) %>% 
  filter(tipoplot == "PRIMARIA") %>% 
  mutate(muestreo = factor(muestreo, levels = c("FEB 20","FEB 21","AGO 21","FEB 22","AGO 22")))
#unique(datos_coords$muestreo)

coords_sf <- leido %>% select(idplot,x,y) %>% filter(!is.na(x)&!is.na(y)) %>% 
  st_as_sf(coords = c("x","y"), crs = 32721) # EPSG 32721 porque en UTM 21S
coords_wgs84 <- st_transform(coords_sf, 4326)


datos_coords_utm <- datos %>%
  left_join(coords_sf,by = "idplot") 
datos_coords_wgs84 <- datos %>% 
  left_join(coords_wgs84,by="idplot")


# ahora sí hago cosas:

#anysign
mapa_any <- ggplot(datos_coords_utm)+
  geom_sf(aes(geometry=geometry,size=anysign_tot))+
  geom_sf(data = datos_coords_utm %>% filter(anysign_tot == 0),
          aes(geometry = geometry), 
          shape = 4, size = 1, stroke = .8, color = "violet") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(~muestreo)

# hozadas
mapa_hoz <- ggplot(datos_coords_utm)+
  geom_sf(aes(geometry=geometry,size=hozadas))+
  geom_sf(data = datos_coords_utm %>% filter(hozadas == 0),
          aes(geometry = geometry), 
          shape = 4, size = 1, stroke = .8, color = "gold") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(~muestreo)

map <- ggplot(datos_coords) +
  geom_sf(aes(color = anysign_tot), size = 3) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Fecha: {frame_time}") +
  transition_time(fecha) +
  ease_aes("linear")


head(datos_coords_utm %>% select("muestreo","idplot","zona",
                                 "geometry","x","y","heces",
                                 "huellas","hozadas","habitat",
                                 "anysign_tot"))

# exploración en el tiempo espacio

resumen <- datos_coords_utm %>%
  group_by(idplot, muestreo) %>%
  summarise(x = first(x),
            y = first(y),
            anysign_tot = max(anysign_tot),
             anysign_tot = max(anysign_tot),
            .groups = "drop") %>%
  arrange(idplot, muestreo) %>%
  group_by(idplot) %>%
  mutate(prev = lag(anysign_tot, default = NA),
         cambio = case_when(
           is.na(prev) ~ "sin_dato",
           anysign_tot > prev ~ "aumentó",
           anysign_tot < prev ~ "disminuyó",
           TRUE ~ "igual"
         ))

# Color personalizado
colores <- c("aumentó" = "red", "disminuyó" = "blue", "igual" = "grey80", "tiempo 0" = "black")

# Figura
ggplot(resumen, aes(x = x, y = y)) +
  geom_point(aes(color = cambio, size = anysign_tot+0.1)) +
  scale_color_manual(values = colores, name = "Cambio") +
  facet_grid(~muestreo) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )+
  labs(title = "Cambios en cantidad de signos entre muestreos",
       x = "", y = "")

# hozadas:
resumen2 <- datos_coords_utm %>%
  group_by(idplot, muestreo) %>%
  summarise(x = first(x),
            y = first(y),
            hozadas = max(hozadas),
            .groups = "drop") %>%
  arrange(idplot, muestreo) %>%
  group_by(idplot) %>%
  mutate(prev = lag(hozadas, default = NA),
         cambio = case_when(
           is.na(prev) ~ "sin_dato",
           hozadas > prev ~ "aumentó",
           hozadas < prev ~ "disminuyó",
           TRUE ~ "igual"
         ))

# Color personalizado
colores <- c("aumentó" = "red", "disminuyó" = "blue", "igual" = "grey80", "tiempo 0" = "black")

# Figura
ggplot(resumen2, aes(x = x, y = y)) +
  geom_point(aes(color = cambio, size = hozadas+0.1)) +
  scale_color_manual(values = colores, name = "Cambio") +
  facet_grid(~muestreo) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )+
  labs(title = "Cambios en hozadas entre muestreos",
       x = "", y = "")


# heatmap tentativo:
ggplot(datos_coords_utm, aes(x = x, y = y)) +
  geom_point(aes(color = anysign_tot,size = hozadas)) +
  scale_color_viridis_c(name = "Signos", option = "C", direction = -1) +
  facet_grid(~ muestreo) +
  coord_equal() +
  theme_minimal()

# intento 2
library(sf)
library(dplyr)
library(spatstat.geom)
library(spatstat.core)
library(terra)
library(ggplot2)
library(purrr)

distmat <- as.matrix(dist(cbind(datos_coords_utm$x,datos_coords_utm$y)))
maxdist <- 2/3*max(distmat) # alrededor de 8500
max(distmat) # 13 mil

# pruebo con conjunto de datos más chico:
# solo 1 muestreo, solo hozadas y xy
idplot <- unique(datos_coords_utm$idplot)
chico <- datos_coords_utm %>% 
  filter(muestreo == "FEB 21") %>% 
  distinct(idplot, .keep_all = TRUE) %>% 
  select(x,y,hozadas)
chico$x <- as.numeric(chico$x)
chico$y <- as.numeric(chico$y)
chico$hozadas <- as.numeric(chico$hozadas)

str(chico)
library(ncf)
correlog <- ncf::correlog(x=chico$x,y=chico$y,z=chico$hozadas, increment = 100,resamp = 99)
plot(correlog)
abline(h=0)

library(geoR)
# fabrioo objeto geoR
objgeoR <- as.geodata(chico)
plot(objgeoR)
emp.variogram <- variog(objgeoR,max.dist = 9000)
plot(emp.variogram)
emp.variogram <- variog4(objgeoR,max.dist = 9000)
plot(emp.variogram)

# kriging
# primero exponential variogram
mlexp <- likfit(objgeoR, cov.model = "exp", ini.cov.pars = c(200,15))
summary(mlexp)
x_vals <- seq(from = min(chico$x), to = max(chico$x), length.out = 25)  # Ajusta a una cantidad razonable
y_vals <- seq(from = min(chico$y), to = max(chico$y), length.out = 25)  # Ajusta a una cantidad razonable
new.grid <- expand.grid(x_vals, y_vals)

#new.grid <- expand.grid(375000:max(chico$x),6467320:max(chico$y))
krig <- krige.conv(objgeoR,locations = new.grid,
                   krige = krige.control(cov.pars = c(mlexp$cov.pars[1],
                                                      mlexp$cov.pars[2]),
                                         nugget=mlexp$nugget,cov.model = "exp",type.krige = "OK"))
image(krig)

library(ggplot2)

library(ggplot2)

# Asegúrate de tener las ubicaciones de la predicción y las predicciones
# krig$prediction.locations es donde están las coordenadas de la cuadrícula
krig_data <- data.frame(
  x = new.grid[, 1],  # Coordenada x
  y = new.grid[, 2],  # Coordenada y
  hozadas = krig$predict               # Las predicciones del kriging
)

# Crear el mapa con ggplot2
ggplot(krig_data, aes(x = x, y = y, fill = hozadas)) +
  geom_tile() +  # Utiliza geom_tile para crear una superficie
  scale_fill_gradient(low = "blue", high = "red") + 
  labs(title = "Kriging: Superficie Interpolada",
       fill = "Hozadas") +  
  theme_minimal()  

