     transfo.fct <- function(Lon,Lat)
      {
tab <- data.frame(Latitude = Lat, Longitude = Lon) 
	#La première étape est de crer un SpatialPointsDataFrame 
        #(équivalent d'un shape sous R) à partir du tableau
	 library(rgdal)
	coords <- as.data.frame(cbind(Lon, Lat)) 
        #petite subtilité, la première colonne contient les Longitudes!
	SP.temp <- SpatialPointsDataFrame(coords = coords ,data=as.data.frame(rep(1,length(Lon))) ,proj4string = CRS("+init=epsg:4326"))
		#comme dans un shape on renseigne dans la fonction 
                #les coordonnées des points (argument coords), 
                #la table attributaire (argument data, non indispensable) 
                #et surtout le sytème de coordonnées (arguement proj4string). 
                #Ici le système géographique est le WGS84, 
                #soit on écrit la formule complète: CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), 
                #soit on écrit le code init: CRS("+init=epsg:4326").
		
	attributes(SP.temp) #l'objet crée est de la classe S4, il s'agit d'un système de listes imbriquées.
	coordinates(SP.temp) 
          #les coordonnées sont toujours en Longitude Latitude, 
          #pour les transformer il faut utiliser la fonction spTransform
	
	SP.proj <-  spTransform(SP.temp, CRS("+init=epsg:3034")) 
            #avec pour arguments le SpatialPointsDataFrame crée précédement et 
            #le système de projection choisi. Dans l'exemple ce doit être de l'ETRS89.
	coordinates(SP.proj) #maintenant les coordonnées sont projetées et tu peux les récupérer.
	
	tab$X <- coordinates(SP.proj)[, 1]
	tab$Y <- coordinates(SP.proj)[, 2]
	return(tab)
}

invtransfo.fct <- function(x,y)
{
 library(rgdal)
# Source data
	coords <- as.data.frame(cbind(x, y)) 
	SP.temp <- SpatialPointsDataFrame(coords =coords ,
     data=as.data.frame(rep(1,length(y))),
     proj4string = CRS("+init=epsg:3034"))
	SP.proj <-  spTransform(SP.temp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
ww <- coordinates(SP.proj)
return(ww)
}

