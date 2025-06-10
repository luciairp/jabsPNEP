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

# Lista con estimaciones kernel por muestreo
densidades_df <- datos_coords_utm %>%
  st_drop_geometry() %>%
  filter(anysign_tot > 0) %>%
  split(.$muestreo) %>%
  map_dfr(function(df) {
    # Crear ventana rectangular
    win <- as.owin(c(range(df$x), range(df$y)))
    pp <- ppp(x = df$x, y = df$y, window = win, marks = df$anysign_tot)
    
    # Expandimos puntos según intensidad
    pp_exp <- pp[rep(1:npoints(pp), marks(pp))]
    
    # Estimar densidad kernel
    dens <- density(pp_exp, sigma = 100) # ajustá sigma si es necesario
    
    # Convertimos a raster y luego a data.frame
    r <- rast(dens)
    df_r <- as.data.frame(r, xy = TRUE)
    names(df_r)[3] <- "densidad"
    df_r$muestreo <- df$muestreo[1]
    df_r
  }, .id = NULL)

# Graficar
ggplot(densidades_df, aes(x = x, y = y, fill = densidad)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Densidad") +
  coord_equal() +
  facet_grid(. ~ muestreo) +
  theme_minimal() +
  labs(title = "Densidad kernel de signos por muestreo",
       x = "X (UTM)", y = "Y (UTM)")
