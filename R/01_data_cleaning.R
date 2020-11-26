#########################
##### DATA CLEANING #####
#########################


### Import the raw data
raw_data <- jk.dusz.tarping::import_raw_data()


### Transform character variables into factors
raw_data <- dplyr::mutate(raw_data,
                          xp_id = as.factor(xp_id),
                          country = as.factor(country),
                          goals = as.factor(goals),
                          planned_duration = factor(x = planned_duration, ordered = TRUE),
                          operation_type = as.factor(operation_type),
                          environment = as.factor(environment),
                          fabric_type = as.factor(fabric_type),
                          season = as.factor(season),
                          preparation = as.factor(preparation),
                          age = factor(x = age, ordered = TRUE),
                          strips_fixation = as.factor(strips_fixation),
                          fabric_fixation = as.factor(fabric_fixation),
                          add_control_type = as.factor(add_control_type),
                          untarped_regrowth = as.factor(untarped_regrowth),
                          latest_condition = as.factor(latest_condition),
                          latest_regrowth = as.factor(latest_regrowth),
                          latest_months = as.factor(latest_months))
# As well as boolean variables:
raw_data <- dplyr::mutate(raw_data,
                          restoration = as.factor(restoration),
                          multiple_ops = as.factor(multiple_ops),
                          liner_geomem = as.factor(liner_geomem),
                          agri_geomem = as.factor(agri_geomem),
                          woven_geotex = as.factor(woven_geotex),
                          mulching_geotex = as.factor(mulching_geotex),
                          pla_geotex = as.factor(pla_geotex),
                          weedsp_geotex = as.factor(weedsp_geotex),
                          other_unknown = as.factor(other_unknown),
                          maxveg = as.factor(maxveg),
                          levelling = as.factor(levelling),
                          fully_tarped = as.factor(fully_tarped),
                          multi_strips = as.factor(multi_strips),
                          tarpfix_multimethod = as.factor(tarpfix_multimethod),
                          trench = as.factor(trench),
                          plantation = factor(x = plantation, ordered = TRUE),
                          repairs = as.factor(repairs),
                          add_control = as.factor(add_control),
                          degradation = as.factor(degradation),
                          pb_fixation = as.factor(pb_fixation),
                          pb_durability = as.factor(pb_durability),
                          pb_trampiercing = as.factor(pb_trampiercing),
                          pb_vandalism = as.factor(pb_vandalism),
                          regrowth_during = as.factor(regrowth_during),
                          reg_staples = as.factor(reg_staples),
                          reg_stripoverlaps = as.factor(reg_stripoverlaps),
                          reg_obstacles = as.factor(reg_obstacles),
                          reg_holes = as.factor(reg_holes),
                          reg_plantations = as.factor(reg_plantations),
                          reg_pierced = as.factor(reg_pierced),
                          reg_edges = as.factor(reg_edges),
                          reg_nearby = as.factor(reg_nearby),
                          tarping_abandoned = as.factor(tarping_abandoned),
                          tarping_completed = as.factor(tarping_completed),
                          tarping_ongoing = as.factor(tarping_ongoing),
                          eff_eradication = as.factor(eff_eradication))
# NOTE: I could probably create a custom function to do this automatically (e.g. for dummy variables,
# with something like: IF variable !is.different(0, 1, NA)... THEN as.factor)? I don't have the time now
# to do that, but it would be a time saver! (I ask if it already exists on StackOverflow).

# IMPORTANT NOTE: the variables "plantation" and "age" are ordinal variables and I have thus coded them
# as such (as an ordered factor). However, it might be preferable, from a statistical point of view,
# to consider it as a numeric variable. The 2nd solution would be more parsimonious (less levels and
# thus lighter models) but would assume that intervals between each level (between 0 and 1,
# between 1 and 2, etc.) are equals when it's not a necessarily true assumption (for "plantation",
# it's probably not). The 1st solution will cause models to add polynomial terms to my factor levels
# (see marked pages in my web browsers).
# I also coded "planned_duration" as an ordinal variable but it there's no problem here because I will
# probably not use it in statistical analyses.


### Description
summary(raw_data)


###


# SHOULD I separate my data into several tibbles (response variables, predictors, others)??? But generally
# or for each model group of models????
