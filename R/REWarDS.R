#' Defined daily dose
#'
#' @description
#' definedDailyDose is a function that uses the World Health Organization's (WHO)
#' defined daily dose method to compute the daily dose and days' supply for prescriptions.
#' This method assumes an average daily consumption of a fixed dose, the defined daily dose
#' (DDD),pre-specified by WHO.
#'
#' @param data           Sample simulated data
#' @param unit           Unit of defined daily dose. Default: 1 unit = 7.5 mg
#' @param dspd_qty       Dispensed quantity: Number of the dispensed tablets to the patient
#'                       at each prescription fill
#' @param strength       Strength of the tablet dispensed in milligrams
#' @param id             Unique patient identification number
#' @param tot_dose_disp  Total dose dispensed:
#'                       dispensed quantity x strength of the tablets dispensed for
#'                       each prescription fill
#' @param Pt_level       When TRUE, the estimated daily dose and days' supply are averaged
#'                       for the patient.
#'
#' @References
#' \insertRef{R}{bibtex}
#' World Health Organization. ATC classification index with DDDs,
#' https://www.whocc.no/atc_ddd_index/?code=B01AA03 (2019, accessed Dec 6th 2020).
#'
#' WHO Collaborating Centre for Drug Statistics Methodology, Guidelines for
#' ATC classification and DDD assignment, 2020. Oslo, 2019.
#'
#' Tanskanen A, Taipale H, Koponen M, et al. Drug exposure in register-based
#' research—An expert-opinion based evaluation of methods. PLoS ONE, 2017; 12: e0184070.
#' DOI: 10.1371/journal.pone.0184070.
#'
#' Sinnott SJ, Polinski JM, Byrne S, et al. Measuring drug exposure: concordance between
#' defined daily dose and days' supply depended on drug class.
#' Journal of Clinical Epidemiology, 2016; 69: 107-113. 2015/07/07.
#' DOI: 10.1016/j.jclinepi.2015.05.026.
#'
#'
#' @return DDD_1_Rx_dose:  Daily dose for prescription
#' @return DDD_1_Rx_DS:    Days' supply for prescription
#' @return DDD_1_Pt_dose:  Average daily dose for patient
#' @return DDD_1_Pt_DS:    Average days' supply for patient
#' @export
#'
#' @details
#' The DDD method can be used for any medication.
#' However, its accuracy has been shown to differ between drug classes.
#'
#' @examples
#' #Patient collects 100 tablets of 5 mg warfarin on January 3rd,
#' #and 100 tablets of 7 mg warfarin on February 1st.
#'
#' ## Using 1 unit DDD:
#' definedDailyDose(data = data, unit = 1, Pt_level = F,id = "ID",
#'                  dspd_qty = "DSPD_QTY", strength = "strength",
#'                  tot_dose_disp = NULL)
#'
#' #tot_dose_disp= 500mg on January 3rd and 700 mg for February 1st.
#' #DDD_1_Rx_dose: 7.5 mg for each prescription fill
#' #DDD_1_Rx_DS is: For Jan 3rd:  500/7.5 = 66.66 day;
#' #                 For Feb 1st: 700/7.5=93.33 days
#'
#' # pt_level can be set as TRUE to get mean values for each patient
#' #DDD_1_Pt_dose: (7.5+ 7.5)/2 = 7.5 mg
#' #DDD_1_Pt_DS: (66.66+93.33)/2 = 79.99 days
#'
#'
definedDailyDose <- function(data, WHO_ddd, dspd_qty, strength, id,
                             tot_dose_disp = NULL, Pt_level = FALSE) {

  if (is.null(tot_dose_disp)) {
    data <-
      data %>%
      mutate(tot_dose_disp = !!sym(dspd_qty) * !!sym(strength))
  }
  else {
    data[ , "tot_dose_disp"] <- data[ , tot_dose_disp]
  }

  data[ , "DDD_Rx_dose"] <- WHO_ddd
  data <-
    data %>%
    mutate(DDD_Rx_DS = tot_dose_disp / (WHO_ddd))
  colnames(data)[colnames(data) == "DDD_Rx_DS"] <- "DDD_Rx_DS"

  if (Pt_level) {
    data <-
      data %>%
      group_by(!!sym(id)) %>%
      mutate(DDD_Pt_dose = mean(WHO_ddd),
             DDD_Pt_DS = mean(tot_dose_disp / (WHO_ddd)))

    colnames(data)[colnames(data) %in% c("DDD_Pt_dose", "DDD_Pt_DS")] <-
      c("DDD_Pt_dose", "DDD_Pt_DS")
  }

  if (! is.null(tot_dose_disp)) data[ , "tot_dose_disp"] <- NULL

  data <- as.data.frame(data)
  return(data)
}




#' Fixed tablet
#'
#' @description
#' fixedTablet is a function that computes the daily dose and days'
#' supply for prescriptions by assuming an average daily consumption of a
#' fixed number of tablets (usually 1) per day by the patient.
#'
#' @param data           Sample simulated data
#' @param tablet         Number of tablets assumed to be consumed by the patient per day.
#'                       Default=1
#' @param dspd_qty       Dispensed quantity: Number of the dispensed tablets to the patient
#'                       at each prescription fill
#' @param strength       Strength of the tablet dispensed in milligrams
#' @param id             Unique patient identification number
#' @param servDate       Date of the prescription fill
#' @param tot_dose_disp  Total dose dispensed:
#'                       dispensed quantity x strength of the tablets dispensed for
#'                       each prescription fill
#' @param Pt_level       When TRUE, the estimated daily dose and days' supply are averaged
#'                       for the patient
#'
#'
#' @References
#' Tanskanen A, Taipale H, Koponen M, et al. Drug exposure in register-based
#' research—An expert-opinion based evaluation of methods. PLoS ONE, 2017; 12: e0184070.
#' DOI: 10.1371/journal.pone.0184070
#'
#' @return tot_dose_disp:            Total dose dispensed at prescription fill:
#'                                   dispensed quantity x strength of the tablet dispensed
#' @return fixed_1_tab_Rx_dose:      Daily dose for prescription
#' @return fixed_1_tab_Rx_DS:        Days' supply for prescription
#' @return fixed_1_tab_Pt_dose:      Average daily dose for patient
#' @return fixed_1_tab_Pt_DS:        Average days' supply for patient
#' @export
#'
#' @details
#' The fixed tablet method can be used for any medication.
#' However, its accuracy has been shown to differ between drug classes.
#'
#' @examples
#' #Patient collects 100 tablets of 5 mg warfarin  on January 3rd,
#' #and 100 tablets of 7 mg warfarin on February 1st.
#'
#' #Assuming consumption of 1 tablet per day:
#' fixedTablet(data = data, tablet = 1, Pt_level = F, id = "ID",
#'             dspd_qty = "DSPD_QTY", strength = "strength",
#'             serv_date = "ServDate", tot_dose_disp = NULL)
#'
#' #tot_dose_disp= 500mg on January 3rd and
#'                 700 mg for February 1st.
#' #fixed_1_tab_Rx_dose: 5 mg for the prescription refill on Jan 3rd, 7 mg for prescription refill on Feb 1st.
#' #fixed_1_tab_Rx_DS is: For Jan 3rd:  500/5= 100 day;  For Feb 1st: 700/7= 100 days
#'
#' #pt_level can be set as TRUE to get mean values for each patient
#' #DDD_1_Pt_dose: (5+ 7)/2 = 6 mg
#' #DDD_1_Pt_DS: (100+100)/2 = 100 days
#'
#'
#'
fixedTablet <- function(data,
                        tablet = 1,
                        dspd_qty,
                        strength,
                        id,
                        serv_date,
                        tot_dose_disp =  NULL,
                        Pt_level = FALSE) {

  if (is.null(tot_dose_disp)) {
    data <-
      data %>%
      mutate(tot_dose_disp = !!sym(dspd_qty) * !!sym(strength))
  }
  else {
    data[ , "tot_dose_disp"] <- data[ , tot_dose_disp]
  }

  data <- data %>%
    group_by(!!sym(id), !!sym(serv_date)) %>%
    mutate(fixed_tab_Rx_dose = sum(tot_dose_disp) / sum(!!sym(dspd_qty)))

  data$fixed_tab_Rx_DS <- as.data.frame(data)[ , dspd_qty]

  if (Pt_level) {
    data <-
      data %>%
      group_by(!!sym(id)) %>%
      mutate(fixed_tab_Pt_dose = mean(fixed_tab_Rx_dose),
             fixed_tab_Pt_DS = mean(fixed_tab_Rx_DS))

    colnames(data)[colnames(data) %in% c("fixed_tab_Pt_dose", "fixed_tab_Pt_DS")] <-
      c(paste0("fixed_", tablet, "_tab_Pt_dose"), paste0("fixed_", tablet, "_tab_Pt_DS"))
  }

  colnames(data)[colnames(data) %in% c("fixed_tab_Rx_dose", "fixed_tab_Rx_DS")] <-
    c(paste0("fixed_", tablet, "_tab_Rx_dose"), paste0("fixed_", tablet, "_tab_Rx_DS"))

  if (! is.null(tot_dose_disp)) data[ , "tot_dose_disp"] <- NULL

  data <- as.data.frame(data)
  return(data)
}






#' Fixed window
#'
#' @description
#' fixedWindow is a function that computes the daily dose and days' supply for
#' prescriptions by assuming a fixed number of days of exposure (usually 90 days)
#' for all patients, reflecting the medication supply policies of most medication insurance
#' plans.
#'
#' @param data           Sample simulated data
#' @param window_length  The number of days that patients' supply of medication is assumed to
#'                       last after each prescription refill. Default= 90 days
#' @param dspd_qty       Dispensed quantity: Number of the dispensed tablets to the patient
#'                       at each prescription fill
#' @param strength       Strength of the tablet dispensed in milligrams
#' @param id             Unique patient identification number
#' @param serv_date      Date of the prescription fill
#' @param tot_dose_disp  Total dose dispensed:
#'                       dispensed quantity x strength of the tablets dispensed for
#'                       each prescription fill
#' @param Pt_level       When TRUE, the estimated daily dose and days' supply are averaged
#'                       for the patient
#'
#' @References
#' Tanskanen A, Taipale H, Koponen M, et al. Drug exposure in register-based
#' research—An expert-opinion based evaluation of methods. PLoS ONE, 2017; 12: e0184070.
#' DOI: 10.1371/journal.pone.0184070.
#'
#' @return tot_dose_disp:                 Total dose dispensed at prescription fill:
#'                                        dispensed quantity x strength of the tablet dispensed
#' @return fixed_window_90_wind_Rx_dose:  Daily dose for prescription
#' @return fixed_90_wind_Rx_DS:           Days' supply for prescription
#' @return fixed_90_wind_Pt_dose:         Average daily dose for patient
#' @return fixed_90_wind_Pt_DS:           Average days' supply for patient
#' @export
#'
#' @examples
#' #Patient collects 100 tablets of 5 mg warfarin  on January 3rd,
#' #and 100 tablets of 7 mg warfarin on February 1st.
#'
#' #Assuming window length of 90 days
#'  fixedWindow(data = data, window_length = 90, id = "ID",
#'              dspd_qty = "DSPD_QTY", strength = "strength",
#'              serv_date = "ServDate", tot_dose_disp =  NULL,
#'              Pt_level = TRUE)
#'
#' #tot_dose_disp = 500mg on January 3rd and
#' #                700 mg for February 1st.
#' #fixed_90_wind_Rx_dose : 500/90 = 5.55 mg  for prescription filled on Jan 3rd;
#'                          700/90=7.77 mg for prescription filled on Feb 1st.
#' #fixed_90_wind_Rx_DS: 90 days for all prescriptions
#'
#' #pt_level can be set as TRUE to get mean values for each patient
#' #fixed_90_wind_Pt_dose : (5.55 + 7.77)/2 = 6.66 mg
#' #fixed_90_wind_Pt_DS: (90 + 90)/2 = 90
#'
#'
#'
fixedWindow <- function(data,
                        window_length = 90,
                        dspd_qty,
                        strength,
                        id,
                        serv_date,
                        tot_dose_disp =  NULL,
                        Pt_level = FALSE) {

  if (is.null(tot_dose_disp)) {
    data <-
      data %>%
      mutate(tot_dose_disp = !!sym(dspd_qty) * !!sym(strength))
  }
  else {
    data[ , "tot_dose_disp"] <- data[ , tot_dose_disp]
  }

  data <- data %>%
    group_by(!!sym(id), !!sym(serv_date)) %>%
    mutate(fixed_wind_Rx_dose = sum(tot_dose_disp) / window_length)

  data$fixed_wind_Rx_DS <- window_length

  if (Pt_level) {
    data <-
      data %>%
      group_by(!!sym(id)) %>%
      mutate(fixed_wind_Pt_dose = mean(fixed_wind_Rx_dose),
             fixed_wind_Pt_DS = mean(fixed_wind_Rx_DS))

    colnames(data)[colnames(data) %in% c("fixed_wind_Pt_dose", "fixed_wind_Pt_DS")] <-
      c(paste0("fixed_", window_length, "_wind_Pt_dose"), paste0("fixed_", window_length, "_wind_Pt_DS"))
  }

  colnames(data)[colnames(data) %in% c("fixed_wind_Rx_dose", "fixed_wind_Rx_DS")] <-
    c(paste0("fixed_", window_length, "_wind_Rx_dose"), paste0("fixed_", window_length, "_wind_Rx_DS"))

  if (! is.null(tot_dose_disp)) data[ , "tot_dose_disp"] <- NULL

  data <- as.data.frame(data)
  return(data)
}








#' REWarDS
#'
#' @description
#' REWarDS (Random Effects Warfarin Days' Supply) is a function that estimates patients'
#' individualized average daily dose and subsequently, days' supply, by fitting a random effects
#' linear regression model to patients' cumulative dose over time.
#' Model parameters include a minimal universally-available set of variables from prescription
#' records.
#'
#' @param data             Sample simulated data
#' @param dspd_qty         Dispensed quantity: Number of tablets dispensed to the patient in
#'                         at each prescription fill
#' @param strength         Strength of the dispensed tablets in milligrams
#' @param id               Unique patient identification number
#' @param gap_handling     Method to handle gaps between prescription fills that are more than
#'                         the permissible gap. Currently, gaps can be handled in three ways:
#'                         1) The “None” method: This is the default and it ignores gaps
#'                         2)The “Initial consecutive Rx” method: Starting from the first
#'                         prescription fill, patients' prescription refills are considered
#'                         until the permissible gap is exceeded. REWarDS uses these initial
#'                         prescription refills to estimates patient's individualized daily
#'                         dose. If the permissible gap is exceeded after the first fill,
#'                         there will only be one prescription for REWaRDS to use, and as of
#'                         now, REWarDS is unable to provide estimates of daily dose based on
#'                         a single prescription.
#'                         3) The “Longest consecutive Rx” method: Looks at all periods with
#'                         consecutive
#'                         prescription refills with gaps between them that do not exceed the
#'                         permissible gap) during the follow up, it then counts the number of
#'                         prescription fills in each period, and picks the period with the
#'                         highest number of prescription fills and estimates the patient's
#'                         average daily dose during that period.
#' @param permissible_gap  Gap (in days) allowed between prescription fills
#' @param serv_date        Date of the prescription fill
#' @param tot_dose_disp    Total dose dispensed:
#'                         dispensed quantity x strength of the tablets dispensed for
#'                         each prescription fill
#' @param Pt_level         When TRUE, the estimated dose and days' supply are averaged
#'                         for the patient
#'
#'
#' @Reference
#' Laird NM and Ware JH. Random-effects models for longitudinal data.
#' Biometrics 1982; 38: 963-974. DOI: 10.2307/2529876
#'
#' @return tot_dose_disp:            Total dose dispensed at prescription fill:
#'                                   dispensed quantity x strength of the tablet dispensed
#' @return REWarDS_avg_daily_dose:   Patient's individualized average daily dose
#' @return REWarDS_Rx_DS:            Days' supply for prescription
#' @return REWarDS_Pt_DS:            Average days' supply for patient
#' @export
#'
#' @details
#' REWarDS has been validated for warfarin. It demonstrated excellent performance that was
#' superior to all current alternative methods for estimating days' supply of warfarin.
#' REWarDS could potentially be used for other medications with variable dosing regimens
#' (e.g. tacrolimus), or in populations with high inter-individual variability in drug
#' clearance (e.g. elderly patients). Validation with cohorts of such patients,
#' or medications other than warfarin, has yet to be done.
#'
#' @examples
#' #Patient collects 100 tablets of 5 mg warfarin  on January 3rd,
#' #and 100 tablets of 7 mg warfarin on February 1st.
#'
#' REWarDS(data = data, id = "ID", dspd_qty = "DSPD_QTY", s
#'        strength = "strength", serv_date = "ServDate",
#'        tot_dose_disp =  NULL, Pt_level = FALSE)
#'
#' #tot_dose_disp= 500mg on January 3rd and
#' #               700 mg for February 1st.
#' #REWarDS_avg_daily_dose: patient's individualized average daily dose obtained from regression analysis
#' #REWarDS_Rx_DS: 500mg/ patient's individualized average daily dose, for Jan 3rd
#' #               700mg/patient's individualized average daily dose , for Feb 1st
#'
#' #Pt_level can be set as TRUE to get mean values for each patient
#' #REWarDS_Pt_DS: average of days' supply on Jan 3rd and Feb 1st
#'
#' #Gap handling method can be specified
#' REWarDS(data = data, id = "ID", dspd_qty = "DSPD_QTY",
#'         strength = "strength", serv_date = "ServDate",
#'         tot_dose_disp =  NULL, Pt_level = TRUE,
#'         gap_handling = "Longest consecutive Rx", permissible_gap = 30)
#'#gap: Gap in number of days between each prescription and the prescription preceding it
#'#Rx_count: Number of prescriptions in each period of consecutive prescriptions until
#'#          the permissible gap is exceeded.

REWarDS <- function(data,
                    dspd_qty,
                    strength,
                    id,
                    gap_handling = "none",  ## "Initial consecutive Rx", "Longest consecutive Rx",
                    permissible_gap = NULL,
                    serv_date,
                    tot_dose_disp =  NULL,
                    Pt_level = FALSE) {

  if (is.null(tot_dose_disp)) {
    data <-
      data %>%
      mutate(tot_dose_disp = !!sym(dspd_qty) * !!sym(strength))
  }
  else {
    data[ , "tot_dose_disp"] <- data[ , tot_dose_disp]
  }

  if (gap_handling == "Initial consecutive Rx") {
    if (is.null(permissible_gap) | permissible_gap <= 0) stop("A positive value is required for permissible_gap in consecutive Rx method.")

    data %>%
      arrange(!!sym(id), !!sym(serv_date)) %>%
      group_by(!!sym(id)) %>%
      mutate(gap = c(0, diff(!!sym(serv_date)))) %>%
      filter(cumsum(gap > permissible_gap) == 0) -> data
  }
  if (gap_handling == "Longest consecutive Rx") {
    if (is.null(permissible_gap) | permissible_gap <= 0) stop("A positive value is required for permissible_gap in longest consecutive Rx method.")

    data %>%
      arrange(!!sym(id), !!sym(serv_date)) %>%
      group_by(!!sym(id)) %>%
      mutate(gap = c(0, diff(!!sym(serv_date))),
             Index = cumsum(gap > permissible_gap)) %>%
      group_by(!!sym(id), Index) %>%
      mutate(Rx_count = n()) %>%
      group_by(!!sym(id)) %>%
      filter(Rx_count == max(Rx_count)) %>%
      select(- c(Index)) -> data
  }


  data_filter <- data %>%
    group_by(!!sym(id)) %>%
    mutate(cum_time = as.numeric(difftime(!!sym(serv_date), (!!sym(serv_date))[1], units = "days")) / 30,
           cum_dose = (cumsum(c(0, tot_dose_disp))[- (n() + 1)]) / 30) %>%
    filter(! n() == 1)
  data_filter <- as.data.frame(data_filter)

  REWarDS_formula <- as.formula(paste0("cum_dose ~ cum_time + (cum_time | ", id, ")"))
  REWarDS_model <- lmer(formula = REWarDS_formula, data = data_filter)

  REWarDS_coef <- coef(REWarDS_model)[[1]]
  REWarDS_coef <- as.data.frame(cbind(id = rownames(REWarDS_coef), REWarDS_coef))
  colnames(REWarDS_coef)[colnames(REWarDS_coef) == "id"] <- id

  data_filter$REWarDS_avg_daily_dose <- REWarDS_coef$cum_time[match(data_filter[ , id], REWarDS_coef[ , id])]
  data_filter$REWarDS_Rx_DS <- data_filter$tot_dose_disp / data_filter$REWarDS_avg_daily_dose

  if (Pt_level) {
    data_filter <-
      data_filter %>%
      group_by(!!sym(id)) %>%
      mutate(REWarDS_Pt_DS = mean(REWarDS_Rx_DS))
  }

  colm_list <- c(id, serv_date, "REWarDS_avg_daily_dose", "REWarDS_Rx_DS")
  if (Pt_level) colm_list <- c(colm_list, "REWarDS_Pt_DS")
  data <- merge(data, data_filter[ , colm_list], all.x = TRUE, by = c(id, serv_date))

  data <- as.data.frame(data)
  return(data)
}


