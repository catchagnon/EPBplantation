#' Estimate Weibull parameters to predict the diameter frequency distribution
#' 
#' This function predicts the parameters of the Weibull distribution used to predict the diameter frequency distribution in unthinned white spruce plantations located in Quebec, Canada
#' based on the percentile-based parameter recovery approach (PCT) as described in Liu et al. 2009
#' (https://doi.org/10.1080/02827580802644599)
#'
#' @param id character. unique id for each plantation
#' @param age numeric. age of plantation (years)
#' @param iqs numeric. station quality index, indicator of the productivity of a given site, based on the
#' mean height of the dominant tree in a a 50-yrs-old plantation. 
#' @param stem.ha numeric. number of stems per hectare
#' @param h numeric. mean height of the dominant tree in the plantation (m)
#' @param mdbh numeric. mean diameter at breast height of trees (cm)
#'
#' @return a dataframe with 5 columns named 'id', location', 'shape', 'scale' and 'stem.ha'
#'
#' corresponding to the id of each plantation, the three parameters used to predict the distribution
#' of the diameters in a white spruce plantation, and the stem density
#'
#' @examples
#' weibull.params.pct(id = "ID", 
#'                    age = 23, 
#'                    iqs = 11.04, 
#'                    stem.ha = 1200, 
#'                    h = 8.29, 
#'                    mdbh = 12)
#'
#' @export
#' @importFrom dplyr %>% 
weibull.params.pct <- function(id, age, iqs, stem.ha, h, mdbh)
{
  stopifnot(is.numeric(age),  age > 0)
  stopifnot(is.numeric(iqs),  iqs > 0)
  stopifnot(is.numeric(stem.ha),  stem.ha > 0)
  stopifnot(is.numeric(h),  h > 0)
  stopifnot(is.numeric(mdbh),  mdbh > 0)
  
  # Calcul le 0eme percentile
  w0 <- 24.291 + (0.105*age) + 0.206*h - (3.311*log(stem.ha))
  # Calcul le le 25eme  percentile
  w25 <- 47.948 - (0.022*age) + (0.291*h) - (0.060*iqs) - (5.007*log(stem.ha))
  # Calcul le 50eme percentile
  w50 <- 52.484 + (0.341*h) - (6.475*log(stem.ha))
  # Calcul du 95eme percentile
  w95 <- 45.185 + (0.185*age) + (0.302*h) + (0.673*iqs) - (5.276*log(stem.ha))
  
  # Calcul et ajout du parametre "location" de la fonction de Weibull
  location <- ((((stem.ha^(1/3))*w0)-w50)/((stem.ha^(1/3))-1))
  # Calcul et ajout du parametre "shape" de la fonction de Weibull
  shape <- 2.343088/((log(w95-location)-(log(w25-location))))
  # Calcul et ajout du parametre "scale" de la fonction de Weibull
  scale <- - (location*gamma(1+1/shape))/gamma(1+2/shape) + sqrt((location/gamma(1+2/shape))^2 * ((gamma(1+1/shape))^2 - gamma(1+2/shape)) + (mdbh^2/gamma(1+2/shape)))
  # Creation d'un dataframe pour accueillir les resultats
  weibull.params <- data.frame("id" = unique(id), "location" = location, "shape" = shape, "scale" = scale, "stem.ha" = stem.ha)
  
  return(weibull.params)
}


#' Estimate diameter frequency distribution
#'
#' This function estimates the diameter frequency distribution in unthinned white spruce plantations using
#' the weibull function parameters 'location', 'scale' and 'shape'.
#'
#' @param dbh numeric. vector of diameters at breast height (dbh) to run the simulation
#' @param weibull.params dataframe. output of function weibull.params.pct()
#'
#' @return a dataframe with number of stem per hectare for each dbh value given as input
#'
#' @examples
#' weibull.params = weibull.params.pct(id = "ID", 
#'                    age = 23, 
#'                    iqs = 11.04, 
#'                    stem.ha = 1200, 
#'                    h = 8.29, 
#'                    mdbh = 12)
#' dhp.dist(dbh = 10:20, weibull.params = weibull.params)
#' 
#' @export
dhp.dist = function(dbh, weibull.params)
{
  stopifnot(is.numeric(dbh), all(dbh > 0))
  stopifnot(all(c("id", "location", "shape", "scale", "stem.ha") %in% names(weibull.params)))
  
  #create a line for each dbh of interest for each plantation
  dbh_for_each_plot <- data.frame("id" = rep(weibull.params$id, each = (max(dbh) - min(dbh) + 1)),
                                  "dbh" = rep(dbh, time = length(weibull.params$id)))
  
  dbh_dist_df <- merge(dbh_for_each_plot, weibull.params)
  
  freq <- with(dbh_dist_df, (shape/scale)*(((dbh-location)/scale)^(shape-1))*exp(-(((dbh-location)/scale)^shape)))
  stem.ha <- freq*dbh_dist_df$stem.ha
  stem.ha[is.nan(stem.ha)] <- 0
  stem.ha <- round(stem.ha, digits = 0)
  
  result_df <- data.frame("ID" = dbh_dist_df$id ,"dbh" = dbh_dist_df$dbh, "stem.ha" = stem.ha)
  #remove 0
  result_df <- result_df %>%
    dplyr::filter(stem.ha > 0)
  return(result_df)
}


#' Estimating tree height 
#' 
#' This function predicts the height of trees in unthinned white spruce plantations based on Auger (2016).
#' Values for parameters SDOM and VPOT must be extracted from https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Auger-Isabelle/Note146.pdf
#' 
#' @param id character. unique id for each plantation
#' @param dbh numeric. diameter at breast height (cm)
#' @param mdbh numeric. mean dbh at breast height of the plantation (cm)
#' @param MAT numeric. mean annual temperature
#' @param SDOM numeric. coefficient associated with bioclimatic subdomain
#' @param VPOT numeric. coefficient associated with potential vegetation
#' @param BA numeric. merchantable basal area (m2/ha)
#'
#' @return a dataframe containing the plantation ID, the dbh distribution and the estimated tree heights
#'
#' @examples
#' height.epb(id = as.character(rep(10:20, each = 10)),
#'            dbh = 10:20, 
#'            mdbh = 12, 
#'            MAT = 2.5, 
#'            SDOM = -0.2749, 
#'            VPOT = -0.1112, 
#'            BA = seq(20, 40, by = 2))
#'   
#' @export 
height.epb = function(id, dbh, mdbh, MAT, SDOM, VPOT, BA)  {
  
  stopifnot(is.numeric(MAT))
  stopifnot(is.numeric(dbh))
  stopifnot(is.numeric(mdbh))
  stopifnot(is.numeric(SDOM))
  stopifnot(is.numeric(VPOT))
  
  height = round((1.3 + (-2.9051 + 0.0252  * BA - 0.2230 *(dbh/mdbh^2) + 
                           0.0789 * MAT + SDOM  + VPOT) * log(dbh + 1) + 
                    (2.2411) * (log(dbh + 1))^2),2)
  
  height_df <- as.data.frame(cbind(id = id,
                     dbh = dbh,
                     height = height))
  return(height_df)
}

#' StatSAW : Modelling lumber product assortment 
#' 
#' Implement StatSAW model as described in Auty et al. 2014 (StatSAW: modelling lumber product assortment using zero-inflated Poisson regression,
#' dx.doi.org/10.1139/cjfr-2013-0500), which was refitted for white and black spruce plantation in 2021.
#'
#' This function predicts the number and volume of lumber products that can be obtained in unthinned white spruce plantation given its dbh (diameter at breast height, cm)
#' frequency distribution. If merchantable volume  (dm3) for each dbh valueis not given, it is calculated using the height and dbh or the trees using the method described in
#' "Équation provenant du mémoire de recherche no 160 (2010): Tarif de cubage, tables de rendement et modèles de croissance pour les plantations d'épinette blanche au Québec
#' https://mffp.gouv.qc.ca/publications/forets/connaissances/recherche/Pregent-Guy/Memoire160.pdf" which is fitted for plantations located in Quebec, Canada.
#'
#' @param id character. unique id for each plantation
#' @param dbh numeric. diameters at breast height (cm) 
#' @param vol.dm3 numeric. (optional) merchantable volume (dm3) for each given dbh value
#' @param stem.ha numeric. number of stem per hectare
#' @param h numeric. (optional) height(m) of trees of each dbh value
#' @param return.volume logical. should the function output include the calculated merchantable volume (dm3)
#'
#' @return a dataframe with number and volume in pmp of lumber products for a single tree of given dbh
#' 
#' The function uses the following six product categories
#' 1. 1x3x08 including 1x3x06, 1x3x07, 1x3x08, 1x3x10, 1x3x12, 1x3x14
#' 2. 1x4x08 including 1x4x06, 1x4x07, 1x4x08, 1x4x10, 1x4x14, 1x4x16, 1x6x08
#' 3. 2x3x08 including 2x3x05, 2x3x06, 2x3x07, 2x3x08
#' 4. 2x4x08 including 2x4x05, 2x4x06, 2x4x07, 2x4x08
#' 5. 2x4x16 including 2x4x10, 2x4x12, 2x4x14, 2x4x16
#' 6. 2x6x10 including 2x6x07, 2x6x08, 2x6x10, 2x6x12, 2x6x14, 2x6x16
#'
#' @examples
#' statSAW.plantation(id = as.character(rep(10:20, each = 10)),
#'                    dbh = 10:20, 
#'                    stem.ha = seq(10,100, by = 10), 
#'                    h = seq(10,15, by = 0.5),
#'                    return.volume = FALSE)
#'
#' @export
statSAW.plantation = function(id, dbh, vol.dm3, stem.ha, h, return.volume = FALSE)
{
  
  stopifnot(is.numeric(dbh), all(dbh > 0))
  stopifnot(is.numeric(stem.ha), all(stem.ha > 0))
  
  if(missing(vol.dm3)) {
    
    if(missing(h)) {
      stop("Please provide a value for h or vol.dm3")
    }
    
    stopifnot(is.numeric(h), all(h > 0))
    
    Ia <- ifelse((9/dbh) > 0.2, 1, 0)
    Ib <- ifelse((9/dbh) > 0.5, 1, 0)
    Ic <- ifelse((9/dbh) > 0.7, 1, 0)
    
    vol.dm3 = ((0.9994 + (0.0159*(9/dbh)) - (0.1071*(9/dbh)^2) - (0.1322*(9/dbh)^3) -
                  (0.5895*Ia*(9/dbh-0.2)^3) - (1.0263*Ib*(9/dbh-0.5)^3) -
                  (1.8282*Ic*(9/dbh-0.7)^3)) * (0.0344*(dbh^1.8329)*(h^1.1793)))
  }
  
  else {
    stopifnot(is.numeric(vol.dm3), all(vol.dm3 > 0))
  }
  
  #prediction categorie 1
  zero_oneby3_08 <- exp(-11.71879)/ (1 + exp(-11.71879))
  count_oneby3_08 <- exp(-3.223201 + (-0.1430564*dbh) + stats::offset(log(vol.dm3)))
  n_oneby3_08 <- (1 - zero_oneby3_08)*count_oneby3_08
  pmp_oneby3_08 <- n_oneby3_08 * 1.96
  
  #prediction categorie 2
  zero_oneby4_08 <- exp(20.350813 + (-1.230441*dbh))/ (1 + exp(20.350813 + (-1.230441*dbh)))
  count_oneby4_08 <- exp(-5.692623449 + (-0.005871257*dbh) + stats::offset(log(vol.dm3)))
  n_oneby4_08 <- (1 - zero_oneby4_08)*count_oneby4_08
  pmp_oneby4_08 <- n_oneby3_08 * 2.87
  
  #prediction categorie 3
  zero_twoby3_08 <- exp(-13.14333)/ (1 + exp(-13.14333))
  count_twoby3_08 <- exp(0.06050313 + (-0.24948850*dbh) + stats::offset(log(vol.dm3)))
  n_twoby3_08 <- (1 - zero_twoby3_08)*count_twoby3_08
  pmp_twoby3_08 <- n_twoby3_08 * 3.49
  
  #prediction categorie 4
  zero_twoby4_08 <- exp(-3.397601 )/ (1 + exp(-3.397601 ))
  count_twoby4_08 <- exp(-3.12273395 + (-0.03829057*dbh) + stats::offset(log(vol.dm3)))
  n_twoby4_08 <- (1 - zero_twoby4_08)*count_twoby4_08
  pmp_twoby4_08 <- n_twoby4_08 * 5.02
  
  #prediction categorie 5
  zero_twoby4_16 <- exp(5.5843517 + (-0.3384668 * dbh))/ (1 + exp(5.5843517 + (-0.3384668 * dbh)))
  count_twoby4_16 <- exp( -3.945173 + stats::offset(log(vol.dm3)))
  n_twoby4_16 <- (1 - zero_twoby4_16)*count_twoby4_16
  pmp_twoby4_16 <- n_twoby4_16 * 9.52
  
  #prediction categorie 6
  zero_twoby6_10 <- exp(22.9506283 + (-0.9861971 * dbh))/ (1 + exp(22.9506283 + (-0.9861971 * dbh)))
  count_twoby6_10 <- exp(-4.9135 + stats::offset(log(vol.dm3)))
  n_twoby6_10 <- (1 - zero_twoby6_10)*count_twoby6_10
  pmp_twoby6_10 <- n_twoby6_10 * 10.5
  
  lumber <- data.frame(id = id, 
                       dbh = dbh,
                       n_oneby3_08 = n_oneby3_08, pmp_oneby3_08 = pmp_oneby3_08,
                       n_oneby4_08 = n_oneby4_08, pmp_oneby4_08 = pmp_oneby4_08,
                       n_twoby3_08 = n_twoby3_08, pmp_twoby3_08 = pmp_twoby3_08,
                       n_twoby4_08 = n_twoby4_08, pmp_twoby4_08 = pmp_twoby4_08,
                       n_twoby4_16 = n_twoby4_16, pmp_twoby4_16 = pmp_twoby4_16,
                       n_twoby6_10 = n_twoby6_10, pmp_twoby6_10 = pmp_twoby6_10)
  
  if (return.volume == TRUE) {
    lumber <- cbind(vol.dm3, lumber)
    return(lumber)
  }
  
  
  else {
    return(lumber)
  }
  
}

#' Estimating the monetary value of lumber product assortments
#' 
#' This function estimate the monetary value of with lumber product assortments
#' @param lumber lumber product assortment. output of function statSAW.plantation().  
#' @param value dataframe comprising two columns; 1. category of lumber product (as in lumber) and 2. associated value
#' 
#' @export
statSAW.value = function(lumber, value) {
  
  pmp_oneby3_08 = NULL
  pmp_oneby4_08 = NULL
  pmp_twoby3_08 = NULL
  pmp_twoby4_08 = NULL
  pmp_twoby4_16 = NULL
  pmp_twoby6_10 = NULL
  id = NULL
  value_tot = NULL
  
  names(value) <- c("product", "price")
  

  lumber = lumber %>% 
      dplyr::group_by(id) %>% 
      dplyr::summarise(id = unique(id),
                pmp_1x3x8 = sum(pmp_oneby3_08),
                pmp_1x4x8 = sum(pmp_oneby4_08),
                pmp_2x3x8 = sum(pmp_twoby3_08),
                pmp_2x4x8 = sum(pmp_twoby4_08),
                pmp_2x4x16 = sum(pmp_twoby4_16),
                pmp_2x6x10 = sum(pmp_twoby6_10))
    
    
    lumber <- reshape2::melt(lumber, id = "id")
    names(lumber) <- c("id", "product", "vol_pmp")
    lumber$product <- substr(lumber$product,5,10)
    
    lumber <- merge(lumber, value)
    lumber$value_tot <- lumber$vol_pmp /1000 * lumber$price
    
   
    lumber <- as.data.frame(lumber %>% 
      dplyr::group_by(id) %>% 
      dplyr::summarise(value = sum(value_tot)))
    
    return(lumber)
  
}         

