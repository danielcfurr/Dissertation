// Change s to run addition sims or in parallel
local s = 1


// Function to simulate data ---------------------------------------------------

capture program drop sim_dataset
program define sim_dataset
	
	syntax [, nitems(integer 32) npersons(integer 500) ///
		tau(real .5) sigma(real 1) seed(integer -1) ]
	
	// Hard coded coefficients
	local beta1 = -.5
	local beta2 = 1
	local beta3 = .5
	local beta4 = .5
	local beta5 = -.5
	
	// Use seed if one is provided
	clear
	if(`seed' > -1) set seed `seed'
	
	// Set up dataset of crossed dummy covariates x3 and x4
	quietly set obs 2
	forvalues i = 3/4 {
		egen x`i' = seq(), from(0) to(1)
	}
	fillin x*
	
	// Expand data to I items and add continuous covariate x2
	generate temp_one = 1
	quietly: expandcl `=`nitems'/_N', cluster(temp_one) generate(temp_item_set)
	sort temp_item_set x3 x4
	generate item = _n
	quietly summarize temp_item_set
	generate x2 = (temp_item_set - 1) / (r(max) - 1) * 2 - 1
	order item x2 x3 x4
	
	// Generate "fixed" and "random" parts of item difficulties
	generate xB = `beta1' + `beta2'*x2 + `beta3'*x3 + `beta4'*x4 + `beta5'*x2*x3
	generate epsilon_training = rnormal(0, `tau')
	generate epsilon_validation = rnormal(0, `tau')
	generate epsilon_test = rnormal(0, `tau')
	
	// Expand dataset to P persons and generate abilities
	quietly: expandcl `npersons', cluster(temp_one) generate(person)
	quietly: generate temp_theta_training = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_training = mean(temp_theta_training)
	quietly: generate temp_theta_validation = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_validation = mean(temp_theta_validation)
	quietly: generate temp_theta_test = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_test = mean(temp_theta_test)
	
	// Simulate responses
	generate y_sameitems = /// "Naive" training data
		rbinomial(1, invlogit(theta_training - xB - epsilon_validation))
	generate y_newitems = /// Better training data
		rbinomial(1, invlogit(theta_training - xB - epsilon_training))
	generate y_validation = ///
		rbinomial(1, invlogit(theta_validation - xB - epsilon_validation))
	generate y_test = ///
		rbinomial(1, invlogit(theta_test - xB - epsilon_test))
	
	sort person item
	drop _fillin temp_*
	
end


// Function to calculate root mean squared error for item predictions ----------

capture program drop rmse
program define rmse, rclass

	syntax ,  [predict_opts(string) predict_var(varname)]
	
	tempvar delta_hat squared_error
	
	if (e(cmd) == "meglm") local if "if person == 1"
	
	// note: delta_hat is easiness, so sign is reversed
	if ("`predict_var'" == "") {
		quietly predict `delta_hat' `if', `predict_opts'
	}
	else {
		quietly generate `delta_hat' = `predict_var' `if'
	}
	
	quietly generate `squared_error' = (-1*`delta_hat' - (xB + epsilon_test))^2
	quietly summarize `squared_error'
	return scalar rmse_full = sqrt(r(mean))
	
	quietly replace `squared_error' = (-1*`delta_hat' - (xB))^2
	quietly summarize `squared_error'
	return scalar rmse_fix = sqrt(r(mean))
	
end


// Function to conduct LLTM analysis -------------------------------------------

capture program drop lltm_cv
program define lltm_cv, rclass
preserve

	syntax varlist (fv), [noquietly]

	if ("`quietly'" == "noquietly") {
		local quietly ""
	}
	else {
		local quietly "quietly"
	}
	
	local converge_fails = 0
	
	// --- Holdout validation ---.
	
	foreach n in newitems sameitems {
		
		`quietly' melogit y_`n' `varlist' || person:
		//rmse, predict_opts(xb)
		//return scalar rmse_full_`n' = r(rmse_full)
		//return scalar rmse_fix_`n' = r(rmse_fix)
		if (e(converged) == 0) local ++converge_fails
		matrix B_`n' = r(table)
		matrix E_`n' = e(b)
		
		// Return parameter estimates
		local cols : colnames B_`n'
		local cols `=subinstr("`cols'", "x", "beta", .)'
		local cols `=subinstr("`cols'", "c.", "", .)'
		local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
		local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
		local cols `=subinstr("`cols'", "_cons", "beta1", 1)'
		local cols `=subinstr("`cols'", "_cons", "sigmasq", 1)'
		forvalues j = 1/`=colsof(B_`n')' {
			local cname = word("`cols'", `j')
			return scalar `=`"`cname'"' + "_est_" + `"`n'"'' = el(B_`n',1,`j')
			return scalar `=`"`cname'"' + "_se_" + `"`n'"'' = el(B_`n',2,`j')
		}
		
		// Get deviance in validation and test data
		local cols : colfullnames E_`n'
		foreach m in validation test {
			matrix E_`n'_`m' = E_`n'
			matrix colnames E_`n'_`m' = `=subinstr("`cols'", "_`n'", "_`m'", .)'
			`quietly' melogit y_`m' `varlist' || person:, asis from(E_`n'_`m') iterate(0)
			return scalar dev_`m'_`n' = -2 * e(ll)
		}
		
	}
	
	// --- Insample deviance (test data only) ---
	
	`quietly' melogit y_test `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_intest = -2 * e(ll)
	
	// --- CV methods (validation data only) ---
	
	`quietly' melogit y_validation `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_invalidation = -2 * e(ll)

	// Information criteria
	`quietly' estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)

	return scalar converge_fails = `converge_fails'
	
restore
end


// Function to conduct two-stage analysis --------------------------------------

capture program drop twostage_cv
program define twostage_cv, rclass
preserve

	syntax varlist (fv), [noquietly]
	
		if ("`quietly'" == "noquietly") {
		local quietly ""
	}
	else {
		local quietly "quietly"
	}
	
	local converge_fails = 0
	`quietly' summarize item
	local I = r(max)

	// Set up meta-analysis dataset
	sort person item
	foreach n in newitems sameitems validation test {
		`quietly' melogit y_`n' ibn.item, noconstant || person:
		if (e(converged) == 0) local ++converge_fails
		if ("`n'" == "validation") return scalar dev_rasch = -2*e(ll)
		matrix A_`n' = r(table)
		matrix B_`n' = A_`n'[1..2, 1..`I']'
		svmat B_`n', names(col)
		rename (b se) (b_`n' se_`n')
	}
	`quietly' drop if b_newitems == .
	
	// Holdout validation meta-regression
	foreach n in newitems sameitems {
		
		`quietly' gsem (b_`n' <- c.se_`n'#c.M@1 `varlist'), variance(M@1)
		matrix B_`n' =  r(table)
		if (e(converged) == 0) local ++converge_fails
		predict delta_hat_`n', eta conditional(fixedonly)
		//rmse, predict_var(delta_hat_`n')
		//return scalar rmse_full_`n' = r(rmse_full)
		//return scalar rmse_fix_`n' = r(rmse_fix)
		
		// Return parameter estimates
		local cols : colnames B_`n'
		local cols `=regexr("`cols'", "c.se_(new|same)items#c.M", "drop")'
		local cols `=subinstr("`cols'", "x", "beta", .)'
		local cols `=subinstr("`cols'", "c.", "", .)'
		local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
		local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
		local cols `=subinstr("`cols'", "_cons", "beta1", 1)'
		local cols `=subinstr("`cols'", "_cons", "drop", 1)'
		local cols `=subinstr("`cols'", "_cons", "tausq", 1)'
		forvalues j = 1/`=colsof(B_`n')' {
			local cname = word("`cols'", `j')
			if ("`cname'" != "drop") {
				return scalar `=`"`cname'"' + "_est_" + `"`n'"'' = el(B_`n',1,`j')
				return scalar `=`"`cname'"' + "_se_" + `"`n'"'' = el(B_`n',2,`j')
			}
		}
		
		matrix E_`n' = e(b)
		generate tausq_`n' = el(E_`n', 1, colsof(E_`n'))
		foreach m in validation test {
			generate ll_`m'_`n' = ///
				log(normalden(b_`m', delta_hat_`n', sqrt(tausq_`n' + se_`m'^2)))
			`quietly' summarize ll_`m'_`n'
			return scalar dev_`m'_`n' = -2*r(sum)
		}
		
	}

	// --- In-sample meta-regression (test data only) ---

	`quietly' gsem (b_test <- c.se_test#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_intest = -2*e(ll)	
	
	// --- CV methods (validation data only) ---
	
	`quietly' gsem (b_validation <- c.se_validation#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_invalidation = -2*e(ll)
	
	// Information criteria
	`quietly' estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)
	
	// LOO CV meta-regression
	quietly generate ll_loo = .
	forvalues i = 1/`I' {
		`quietly' gsem (b_validation <- c.se_validation#c.M@1 `varlist') ///
			if item != `i', variance(M@1)
		if (e(converged) == 0) local ++converge_fails
		predict delta_hat_temp, eta conditional(fixedonly)
		matrix E_loo = e(b)
		local tau2 = el(E_loo, 1, colsof(E_loo))
		`quietly' replace ll_loo = log(normalden(delta_hat_temp, b_validation, ///
			sqrt(`tau2' + se_validation^2))) if item == `i'
		drop delta_hat_temp
	}
	`quietly' summarize ll_loo
	return scalar dev_loo = -2*r(sum)
	
	// TRMSEA
	//generate ll_saturated = log(normalden(0, 0, se))
	//quietly summarize ll_saturated
	//local dev_saturated = -2*r(sum)
	//local df = e(N) - e(rank)
	//return scalar rmsea = ///
	//	sqrt(1/(e(N)*`df') * max(`dev_insample' - `dev_saturated' - `df', 0))
	
	return scalar converge_fails = `converge_fails'
	
restore
end


// Simulation ------------------------------------------------------------------

// Import random number seeds
import delimited "random.txt", delimiter(tab) varnames(nonames) clear
rename v`s' seed
drop v*
keep in 1/500

// Set up simulation conditions, store in a matrix
egen condition = seq(), from(1) to(5)
generate double tau = .5
replace tau = .2 if condition == 1
replace tau = 1 if condition == 3
generate nitems = 32
replace nitems = 16 if condition == 4
replace nitems = 64 if condition == 5
generate npersons = 500
mkmat *, matrix(SIM)
clear

// Define models (predictors)
local m1 "x2-x4"
local m2 "x2-x4 c.x2#c.x3"
local m3 "x2-x4 c.x2#c.x3 c.x2#c.x4"

// Prepare to store results
postfile memhold_lltm double(    ///
		seed                     ///
		condition                ///
		tau                      ///
		nitems                   ///
		npersons                 ///
		model                    ///
		dev_invalidation         ///
		aic                      ///
		bic                      ///
		dev_intest               ///
		dev_test_newitems        ///
		dev_validation_newitems  ///
		dev_test_sameitems       ///
		dev_validation_sameitems ///
		beta1_est_newitems       ///
		beta2_est_newitems       ///
		beta3_est_newitems       ///
		beta4_est_newitems       ///
		beta5_est_newitems       ///
		beta6_est_newitems       ///
		sigmasq_est_newitems       ///
		beta1_est_sameitems      ///
		beta2_est_sameitems      ///
		beta3_est_sameitems      ///
		beta4_est_sameitems      ///
		beta5_est_sameitems      ///
		beta6_est_sameitems      ///
		sigmasq_est_sameitems      ///
		beta1_se_newitems        ///
		beta2_se_newitems        ///
		beta3_se_newitems        ///
		beta4_se_newitems        ///
		beta5_se_newitems        ///
		beta6_se_newitems        ///
		sigmasq_se_newitems        ///
		beta1_se_sameitems       ///
		beta2_se_sameitems       ///
		beta3_se_sameitems       ///
		beta4_se_sameitems       ///
		beta5_se_sameitems       ///
		beta6_se_sameitems       ///
		sigmasq_se_sameitems       ///
	) using results_lltm_`s', replace
postfile memhold_twostage double(    ///
		seed                     ///
		condition                ///
		tau                      ///
		nitems                   ///
		npersons                 ///
		model                    ///
		dev_invalidation         ///
		aic                      ///
		bic                      ///
		dev_intest               ///
		dev_test_newitems        ///
		dev_validation_newitems  ///
		dev_test_sameitems       ///
		dev_validation_sameitems ///
		beta1_est_newitems       ///
		beta2_est_newitems       ///
		beta3_est_newitems       ///
		beta4_est_newitems       ///
		beta5_est_newitems       ///
		beta6_est_newitems       ///
		tausq_est_newitems       ///
		beta1_est_sameitems      ///
		beta2_est_sameitems      ///
		beta3_est_sameitems      ///
		beta4_est_sameitems      ///
		beta5_est_sameitems      ///
		beta6_est_sameitems      ///
		tausq_est_sameitems      ///
		beta1_se_newitems        ///
		beta2_se_newitems        ///
		beta3_se_newitems        ///
		beta4_se_newitems        ///
		beta5_se_newitems        ///
		beta6_se_newitems        ///
		tausq_se_newitems        ///
		beta1_se_sameitems       ///
		beta2_se_sameitems       ///
		beta3_se_sameitems       ///
		beta4_se_sameitems       ///
		beta5_se_sameitems       ///
		beta6_se_sameitems       ///
		tausq_se_sameitems       ///
		dev_loo                  ///
		dev_rasch                ///
	) using results_twostage_`s', replace
postfile memhold_errors double(seed) condition double(tau) nitems npersons ///
	model str12(function) using results_errors_`s', replace

timer clear
timer on 1
forvalues i = 1/`=rowsof(SIM)' {

	display "$S_TIME: Starting `i' of `=rowsof(SIM)'"

	sim_dataset, npersons(500) sigma(1) ///
		seed(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
		nitems(`=el(SIM, `i', colnumb(SIM, "nitems"))') /// 
		tau(`=el(SIM, `i', colnumb(SIM, "tau"))')
		
	forvalues m = 1/3 {
	
		capture noisily lltm_cv `m`m''
		// If the function fails, either a positive number or nothing will be
		// posted to r(converge_fails). If so, post conditions to 
		// memhold_errors. Otherwise, post simulation results to memhold_lltm.
		if r(converge_fails) != 0 {
			post memhold_errors ///
				(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
				(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
				(`=el(SIM, `i', colnumb(SIM, "tau"))') ///
				(`=el(SIM, `i', colnumb(SIM, "nitems"))') ///
				(`=el(SIM, `i', colnumb(SIM, "npersons"))') ///
				(`m') ///
				("lltm_cv")
		} 
		else {
			post memhold_lltm ///
				(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
				(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
				(`=el(SIM, `i', colnumb(SIM, "tau"))') ///
				(`=el(SIM, `i', colnumb(SIM, "nitems"))') ///
				(`=el(SIM, `i', colnumb(SIM, "npersons"))') ///
				(`m') ///
				(r(dev_invalidation         )) ///
				(r(aic                      )) ///
				(r(bic                      )) ///
				(r(dev_intest               )) ///
				(r(dev_test_newitems        )) ///
				(r(dev_validation_newitems  )) ///
				(r(dev_test_sameitems       )) ///
				(r(dev_validation_sameitems )) ///
				(r(beta1_est_newitems           )) ///
				(r(beta2_est_newitems           )) ///
				(r(beta3_est_newitems           )) ///
				(r(beta4_est_newitems           )) ///
				(r(beta5_est_newitems           )) ///
				(r(beta6_est_newitems           )) ///
				(r(sigmasq_est_newitems         )) ///
				(r(beta1_est_sameitems          )) ///
				(r(beta2_est_sameitems          )) ///
				(r(beta3_est_sameitems          )) ///
				(r(beta4_est_sameitems          )) ///
				(r(beta5_est_sameitems          )) ///
				(r(beta6_est_sameitems          )) ///
				(r(sigmasq_est_sameitems        )) ///
				(r(beta1_se_newitems           )) ///
				(r(beta2_se_newitems           )) ///
				(r(beta3_se_newitems           )) ///
				(r(beta4_se_newitems           )) ///
				(r(beta5_se_newitems           )) ///
				(r(beta6_se_newitems           )) ///
				(r(sigmasq_se_newitems         )) ///
				(r(beta1_se_sameitems          )) ///
				(r(beta2_se_sameitems          )) ///
				(r(beta3_se_sameitems          )) ///
				(r(beta4_se_sameitems          )) ///
				(r(beta5_se_sameitems          )) ///
				(r(beta6_se_sameitems          )) ///
				(r(sigmasq_se_sameitems        ))
		}

		capture noisily twostage_cv `m`m''
		// If the function fails, either a positive number or nothing will be
		// posted to r(converge_fails). If so, post conditions to 
		// memhold_errors. Otherwise, post simulation results to 
		// memhold_twostage.
		if r(converge_fails) != 0 {
			post memhold_errors ///
				(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
				(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
				(`=el(SIM, `i', colnumb(SIM, "tau"))') ///
				(`=el(SIM, `i', colnumb(SIM, "nitems"))') ///
				(`=el(SIM, `i', colnumb(SIM, "npersons"))') ///
				(`m') ///
				("twostage_cv")
		} 
		else {
			post memhold_twostage ///
				(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
				(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
				(`=el(SIM, `i', colnumb(SIM, "tau"))') ///
				(`=el(SIM, `i', colnumb(SIM, "nitems"))') ///
				(`=el(SIM, `i', colnumb(SIM, "npersons"))') ///
				(`m') ///
				(r(dev_invalidation         )) ///
				(r(aic                      )) ///
				(r(bic                      )) ///
				(r(dev_intest               )) ///
				(r(dev_test_newitems        )) ///
				(r(dev_validation_newitems  )) ///
				(r(dev_test_sameitems       )) ///
				(r(dev_validation_sameitems )) ///
				(r(beta1_est_newitems           )) ///
				(r(beta2_est_newitems           )) ///
				(r(beta3_est_newitems           )) ///
				(r(beta4_est_newitems           )) ///
				(r(beta5_est_newitems           )) ///
				(r(beta6_est_newitems           )) ///
				(r(tausq_est_newitems           )) ///
				(r(beta1_est_sameitems          )) ///
				(r(beta2_est_sameitems          )) ///
				(r(beta3_est_sameitems          )) ///
				(r(beta4_est_sameitems          )) ///
				(r(beta5_est_sameitems          )) ///
				(r(beta6_est_sameitems          )) ///
				(r(tausq_est_sameitems          )) ///
				(r(beta1_se_newitems           )) ///
				(r(beta2_se_newitems           )) ///
				(r(beta3_se_newitems           )) ///
				(r(beta4_se_newitems           )) ///
				(r(beta5_se_newitems           )) ///
				(r(beta6_se_newitems           )) ///
				(r(tausq_se_newitems           )) ///
				(r(beta1_se_sameitems          )) ///
				(r(beta2_se_sameitems          )) ///
				(r(beta3_se_sameitems          )) ///
				(r(beta4_se_sameitems          )) ///
				(r(beta5_se_sameitems          )) ///
				(r(beta6_se_sameitems          )) ///
				(r(tausq_se_sameitems          )) ///
				(r(dev_loo))                   ///
				(r(dev_rasch))
		}
	}

}
timer off 1
timer list

postclose memhold_twostage
postclose memhold_lltm
postclose memhold_errors


// Save an example dataset
clear
sim_dataset
save example, replace
