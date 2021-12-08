#' Estimate Weibull parameters to predict the diameter frequency distribution
#' 
#' This function predicts the parameters of the Weibull distribution used to predict the diameter frequency distribution in unthinned white spruce plantations located in Quebec, Canada
#' based on the percentile-based parameter recovery approach (PCT) as described in Liu et al. 2009
#' (https://doi.org/10.1080/02827580802644599)
#'
#' @param age numeric. Age of plantation
#' @param iqs numeric. IQS (Station quality index) indicator of the productivity of a given site, based on the
#' mean height of the dominant tree in a a 50-yrs-old plantation.
#' @param stem_ha numeric. Number of stem per hectare
#' @param h numeric. Mean height of the dominant tree in the plantation
#' @param id character. Unique id for each plantation
#' @param mdbh numeric. Arithmetic mean diameter at breast height of trees.
#'
#' @return a dataframe with 5 columns named 'id', location', 'shape', 'scale' and 'stem_ha'
#'
#' corresponding to the id of each plantation, the three parameters used to predict the distribution
#' of the diameters in a white spruce plantation, and the stem density
#'
#' @examples
#' weibull_param_pct(age = 23, iqs = 11.04, tige_ha = 1200, h = 8.29)
#'
#' @export
weibull.param.pct <- function(age, iqs, stem_ha, h, mdbh, id)
{
  stopifnot(is.numeric(age),  age > 0)
  stopifnot(is.numeric(iqs),  iqs > 0)
  stopifnot(is.numeric(stem_ha),  stem_ha > 0)
  stopifnot(is.numeric(h),  h > 0)
  stopifnot(is.numeric(mdbh),  mdbh > 0)
  
  # Calcul le 0eme percentile
  w0 <- 24.291 + (0.105*age) + 0.206*h - (3.311*log(stem_ha))
  # Calcul le le 25eme  percentile
  w25 <- 47.948 - (0.022*age) + (0.291*h) - (0.060*iqs) - (5.007*log(stem_ha))
  # Calcul le 50eme percentile
  w50 <- 52.484 + (0.341*h) - (6.475*log(stem_ha))
  # Calcul du 95eme percentile
  w95 <- 45.185 + (0.185*age) + (0.302*h) + (0.673*iqs) - (5.276*log(stem_ha))
  
  # Calcul et ajout du parametre "location" de la fonction de Weibull
  location <- ((((stem_ha^(1/3))*w0)-w50)/((stem_ha^(1/3))-1))
  # Calcul et ajout du parametre "shape" de la fonction de Weibull
  shape <- 2.343088/((log(w95-location)-(log(w25-location))))
  # Calcul et ajout du parametre "scale" de la fonction de Weibull
  scale <- - (location*gamma(1+1/shape))/gamma(1+2/shape) + sqrt((location/gamma(1+2/shape))^2 * ((gamma(1+1/shape))^2 - gamma(1+2/shape)) + (mdbh^2/gamma(1+2/shape)))
  # Creation d'un dataframe pour accueillir les resultats
  weibull_params <- data.frame("id" = id, "location" = location, "shape" = shape, "scale" = scale, "stem_ha" = stem_ha)
  
  return(weibull_params)
}


#' Estimate diameter frequency distribution
#'
#' This function estimates the diameter frequency distribution in unthinned white spruce plantations using
#' the weibull function parameters 'location', 'scale' and 'shape'.
#'
#' @param dbh numeric. vector of diameters at breast height (dbh) to run the simulation
#' @param weibull_param dataframe. output of function param_weibull_pct
#'
#' @return a dataframe with number of stem per hectare for each dbh value given as input
#'
#'
#' @export
dhp.dist = function(dbh, weibull_params)
{
  stopifnot(is.numeric(dbh), all(dbh > 0))
  stopifnot(all(c("id", "location", "shape", "scale", "stem_ha") %in% names(weibull_params)))
  
  #create a line for each dbh of interest for each plantation
  dbh_for_each_plot <- data.frame("id" = rep(weibull_params$id, each = (max(dbh) - min(dbh) + 1)),
                                  "dbh" = rep(dbh, time = length(id)))
  
  dbh_dist_df <- merge(dbh_for_each_plot, weibull_params)
  
  freq <- with(dbh_dist_df, (shape/scale)*(((dbh-location)/scale)^(shape-1))*exp(-(((dbh-location)/scale)^shape)))
  stem_ha <- freq*dbh_dist_df$stem_ha
  stem_ha[is.nan(stem_ha)] <- 0
  stem_ha <- round(stem_ha, digits = 0)
  
  result_df <- data.frame("ID" = dbh_dist_df$id ,"dbh" = dbh_dist_df$dbh, "stem_ha" = stem_ha)
  #remove 0
  result_df <- result_df %>%
    filter(stem_ha > 0)
  return(result_df)
}




#' StatSAW : Modelling lumber product assortment 
#' 
#' Implement StatSAW model as described in Auty et al. 2014 (StatSAW: modelling lumber product assortment using zero-inflated Poisson regression,
#' dx.doi.org/10.1139/cjfr-2013-0500), which was refitted for white and black spruce plantation in 2021.
#'
#' This function predicts the number and volume of lumber products that can be obtained in unthinned white spruce plantation given its dbh (diameter at breast height)
#' frequency distribution. If merchantable volume  (dm3) for each dbh valueis not given, it is calculated using the height and dbh or the trees using the method described in
#' "Équation provenant du mémoire de recherche no 160 (2010): Tarif de cubage, tables de rendement et modèles de croissance pour les plantations d'épinette blanche au Québec
#' https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Pregent-Guy/Memoire160.pdf" which is fitted for plantations located in Quebec, Canada.
#'
#'
#' @param dbh numeric. diameters at breast height (dbh) 
#' @param vol_m numeric. (optional) merchantable volume (dm3) for each given dbh value
#' @param stem_ha numeric. Number of stem per hectare
#' @param h numeric. (optional) Height(m) of trees of each dbh value
#'
#' @return a dataframe with number and volume in pmp of lumber products for a diameter frequency distribution
#' 
#' The function uses the following six product categories
#' 1x3x08 <- c("1x3x06", "1x3x07", "1x3x08", "1x3x10", "1x3x12", "1x3x14")
#' 1x4x08 <- c("1x4x06", "1x4x07", "1x4x08", "1x4x10", "1x4x14", "1x4x16", "1x6x08")
#' 2x3x08 <- c("2x3x05", "2x3x06", "2x3x07", "2x3x08")
#' 2x4x08 <- c("2x4x05", "2x4x06", "2x4x07", "2x4x08")
#' 2x4x16 <- c("2x4x10", "2x4x12", "2x4x14", "2x4x16")
#' 2x6x10 <- c("2x6x07", "2x6x08", "2x6x10", "2x6x12", "2x6x14", "2x6x16")
#'
#' @examples
#' param = param_weibull_pct(age = 23, iqs = 11.04, tige_ha = 1200, h = 8.29)
#' stem_per_ha = dhp_dist(1:60, param)
#'
#' @export
statSAW.plantation = function(dbh, vol_dm3, stem_ha, h, return.volume = FALSE)
{
  stopifnot(is.numeric(dbh), all(dbh > 0))
  stopifnot(is.numeric(stem_ha), all(stem_ha > 0))
  
  if(missing(vol_dm3)) {
    
    if(missing(h)) {
      stop("Please provide a value for h or vol_dm3")
    }
    
    stopifnot(is.numeric(h), all(h > 0))
    
    Ia <- ifelse((9/dbh) > 0.2, 1, 0)
    Ib <- ifelse((9/dbh) > 0.5, 1, 0)
    Ic <- ifelse((9/dbh) > 0.7, 1, 0)
    
    vol_dm3 = ((0.9994 + (0.0159*(9/dbh)) - (0.1071*(9/dbh)^2) - (0.1322*(9/dbh)^3) -
                  (0.5895*Ia*(9/dbh-0.2)^3) - (1.0263*Ib*(9/dbh-0.5)^3) -
                  (1.8282*Ic*(9/dbh-0.7)^3)) * (0.0344*(dbh^1.8329)*(h^1.1793)))
  }
  
  else {
    stopifnot(is.numeric(vol_dm3), all(vol_dm3 > 0))
  }
  
  #prediction categorie 1
  zero_oneby3_08 <- exp(-11.71879)/ (1 + exp(-11.71879))
  count_oneby3_08 <- exp(-3.223201 + (-0.1430564*dbh) + offset(log(vol_dm3)))
  n_oneby3_08 <- (1 - zero_oneby3_08)*count_oneby3_08
  pmp_oneby3_08 <- n_oneby3_08 * 1.96
  
  #prediction categorie 2
  zero_oneby4_08 <- exp(20.350813 + (-1.230441*dbh))/ (1 + exp(20.350813 + (-1.230441*dbh)))
  count_oneby4_08 <- exp(-5.692623449 + (-0.005871257*dbh) + offset(log(vol_dm3)))
  n_oneby4_08 <- (1 - zero_oneby4_08)*count_oneby4_08
  pmp_oneby4_08 <- n_oneby3_08 * 2.87
  
  #prediction categorie 3
  zero_twoby3_08 <- exp(-13.14333)/ (1 + exp(-13.14333))
  count_twoby3_08 <- exp(0.06050313 + (-0.24948850*dbh) + offset(log(vol_dm3)))
  n_twoby3_08 <- (1 - zero_twoby3_08)*count_twoby3_08
  pmp_twoby3_08 <- n_twoby3_08 * 3.49
  
  #prediction categorie 4
  zero_twoby4_08 <- exp(-3.397601 )/ (1 + exp(-3.397601 ))
  count_twoby4_08 <- exp(-3.12273395 + (-0.03829057*dbh) + offset(log(vol_dm3)))
  n_twoby4_08 <- (1 - zero_twoby4_08)*count_twoby4_08
  pmp_twoby4_08 <- n_twoby4_08 * 5.02
  
  #prediction categorie 5
  zero_twoby4_16 <- exp(5.5843517 + (-0.3384668 * dbh))/ (1 + exp(5.5843517 + (-0.3384668 * dbh)))
  count_twoby4_16 <- exp( -3.945173 + offset(log(vol_dm3)))
  n_twoby4_16 <- (1 - zero_twoby4_16)*count_twoby4_16
  pmp_twoby4_16 <- n_twoby4_16 * 9.52
  
  #prediction categorie 6
  zero_twoby6_10 <- exp(22.9506283 + (-0.9861971 * dbh))/ (1 + exp(22.9506283 + (-0.9861971 * dbh)))
  count_twoby6_10 <- exp(-4.9135 + offset(log(vol_dm3)))
  n_twoby6_10 <- (1 - zero_twoby6_10)*count_twoby6_10
  pmp_twoby6_10 <- n_twoby6_10 * 10.5
  
  lumber <- data.frame(dbh = dbh,
                       n_oneby3_08 = n_oneby3_08*stem_ha, pmp_oneby3_08 = pmp_oneby3_08*stem_ha,
                       n_oneby4_08 = n_oneby4_08*stem_ha, pmp_oneby4_08 = pmp_oneby4_08*stem_ha,
                       n_twoby3_08 = n_twoby3_08*stem_ha, pmp_twoby3_08 = pmp_twoby3_08*stem_ha,
                       n_twoby4_08 = n_twoby4_08*stem_ha, pmp_twoby4_08 = pmp_twoby4_08*stem_ha,
                       n_twoby4_16 = n_twoby4_16*stem_ha, pmp_twoby4_16 = pmp_twoby4_16*stem_ha,
                       n_twoby6_10 = n_twoby6_10*stem_ha, pmp_twoby6_10 = pmp_twoby6_10*stem_ha)
  
  if (return.volume == TRUE) {
    lumber <- cbind(vol_dm3, lumber)
    return(lumber)
  }
  
  
  else {
    return(lumber)
  }
  
}

#' Estimating the monetary value of lumber product assortments
#' 
#' This function estimate the monetary value of with lumber product assortments
#' @param lumber lumber product assortment. output of function statSAWplantation().  
#' @param value dataframe comprising two columns; 1. category of lumber product (as in lumber) and 2. associated value
#' @param id (optional) Unique id for each plantation, must be the same length as lumber. If missing, a single moneraty value will be calculated for the whole lumber product assortment
#' @export
statSAW.value = function(lumber, value, id = NULL) {
  
  stopifnot(length(id) == length(lumber[,1]))
  names(value) <- c("product", "price")
  
  if (is.null(id)) {
    lumber = lumber %>% 
      summarise(pmp_1x3x8 = sum(pmp_oneby3_08),
                pmp_1x4x8 = sum(pmp_oneby4_08),
                pmp_2x3x8 = sum(pmp_twoby3_08),
                pmp_2x4x8 = sum(pmp_twoby4_08),
                pmp_2x4x16 = sum(pmp_twoby4_16),
                pmp_2x6x10 = sum(pmp_twoby6_10))
    
    lumber <- data.frame(vol_pmp = t(lumber),
                         product = substr(colnames(lumber),5,10))
    
    stopifnot(all(unique(lumber$product) %in% value$product))
    
    lumber <- merge(lumber, value)
    return(sum(lumber$vol_pmp /1000 * lumber$price))
  }
  
  else {
    lumber <- cbind(id, lumber)
    lumber = lumber %>% 
      group_by(id) %>% 
      summarise(id = unique(id),
                pmp_1x3x8 = sum(pmp_oneby3_08),
                pmp_1x4x8 = sum(pmp_oneby4_08),
                pmp_2x3x8 = sum(pmp_twoby3_08),
                pmp_2x4x8 = sum(pmp_twoby4_08),
                pmp_2x4x16 = sum(pmp_twoby4_16),
                pmp_2x6x10 = sum(pmp_twoby6_10))
    
    lumber <- data.frame(vol_pmp = t(lumber),
                         product = substr(colnames(lumber),5,10))
    lumber <- merge(lumber, value)
    lumber <- lumber %>% 
      group_by(id) %>% 
      summarise(value = sum(lumber$vol_pmp /1000 * lumber$price))
    
    return(lumber)
  }
  
}         


#' Estimating tree height 
#' 
#' This function predicts the height of trees in unthinned white spruce plantations based on Auger 2016
#'  (Une nouvelle relation hauteur-diamètre tenant compte de l’influence de la station et du climat pour 27 essences commerciales du Québec) 
#'  https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Auger-Isabelle/Note146.pdf
#' 
#' @param dbh numeric. 
#' @param mdbh numeric. 
#' @param MAT numeric. Mean annual temperature
#' @param SDOM numeric. Coefficient associated with bioclimatic subdomain. See https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Auger-Isabelle/Note146.pdf
#' @param VPOT numeric. Coefficient associated with potential vegetation. See https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Auger-Isabelle/Note146.pdf
#' @param BA numeric. Merchantable basal area (m2/ha)
#'
#' @return a dataframe containing the height of each given dbh value
#'
#' @export 
height.epb = function(dbh, mdbh, MAT, SDOM, VPOT, BA)  {
  stopifnot(is.numeric(MAT))
  stopifnot(is.numeric(dbh))
  stopifnot(is.numeric(mdbh))
  stopifnot(is.numeric(SDOM))
  stopifnot(is.numeric(VPOT))
  
  height = round((1.3 + (-2.9051 + 0.0252  * BA - 0.2230 *(dbh/mdbh^2) + 
                           0.0789 * MAT + SDOM  + VPOT) * log(dbh + 1) + 
                    (2.2411) * (log(dbh + 1))^2),2)
  
  return(height)
}
