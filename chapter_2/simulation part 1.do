// Change s to run addition sims or in parallel
local s = 1


// Function to simulate data ---------------------------------------------------

capture program drop sim_dataset
program define sim_dataset, rclass
	
	syntax [, nitems(integer 32) npersons(integer 500) ///
		 sigma(real 1.5) tau(real .5) b(real .5) seed(integer -1) ]
	
	// Coefficients
	local beta1 = 0
	local beta2 = `b'
	local beta3 = `b'
	local beta4 = `b'
	
	// Use seed if one is provided
	clear
	if(`seed' > -1) set seed `seed'
	
	// Set up dataset of item covariates
	quietly set obs `nitems'
	generate item = _n
	forvalues i = 2/5 {
		generate x`i' = rnormal()
	}
	
	// Generate "fixed" and "random" parts of item difficulties
	generate xB = `beta1' + `beta2'*x2 + `beta3'*x3 + `beta4'*x4
	generate epsilon_training = rnormal(0, `tau')
	generate epsilon_validation = rnormal(0, `tau')
	generate epsilon_test = rnormal(0, `tau')
	
	// Expand dataset to P persons and generate abilities
	generate temp_one = 1
	quietly: expandcl `npersons', cluster(temp_one) generate(person)
	quietly: generate temp_theta_training = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_training = mean(temp_theta_training)
	quietly: generate temp_theta_validation = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_validation = mean(temp_theta_validation)
	quietly: generate temp_theta_test = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta_test = mean(temp_theta_test)
	
	// Simulate responses
	generate y_train = /// Training response data
		rbinomial(1, invlogit(theta_training - xB - epsilon_training))
	generate y_newitems = /// Validation data with new items
		rbinomial(1, invlogit(theta_validation - xB - epsilon_validation))
	generate y_sameitems = /// Validation data with same items
		rbinomial(1, invlogit(theta_validation - xB - epsilon_training))
	generate y_test = /// Test data
		rbinomial(1, invlogit(theta_test - xB - epsilon_test))	
	
	// Set up meta-analysis dataset
	local converge_fails = 0
	quietly summarize item
	local I = r(max)
	sort person item
	foreach n in train newitems sameitems test {
		quietly melogit y_`n' ibn.item, noconstant || person:
		if (e(converged) == 0) local ++converge_fails
		// Return training sample Rasch deviance
		if ("`n'" == "train") return scalar dev_rasch = -2*e(ll)
		matrix A_`n' = r(table)
		matrix B_`n' = A_`n'[1..2, 1..`I']'
		svmat B_`n', names(col)
		rename (b se) (b_`n' se_`n')
	}
	
	// Return training sample "Rasch" meta-regression deviance
	quietly generate ll_train = log(normalden(0, 0, se_train))
	quietly summarize ll_train
	return scalar dev_meta = -2 * r(sum)
			
	return scalar converge_fails = `converge_fails'
	
	sort person item
	drop temp_*
	
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
	
	// Fit to training data
	`quietly' melogit y_train `varlist' || person:
	return scalar dev_train = -2 * e(ll)
	if (e(converged) == 0) local ++converge_fails
	matrix B = r(table)
	matrix E = e(b)
	
	// Get RMSE for HV
	rmse, predict_opts(xb)
	return scalar rmse_full = r(rmse_full)
	return scalar rmse_fix = r(rmse_fix)
	
	// Information criteria
	`quietly' estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)
	
	// Return parameter estimates
	local cols : colnames B
	local cols `=subinstr("`cols'", "x", "beta", .)'
	local cols `=subinstr("`cols'", "c.", "", .)'
	//local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
	//local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
	local cols `=subinstr("`cols'", "_cons", "beta1", 1)'
	local cols `=subinstr("`cols'", "_cons", "sigmasq", 1)'
	forvalues j = 1/`=colsof(B)' {
		local cname = word("`cols'", `j')
		return scalar `cname'_est = el(B,1,`j')
		return scalar `cname'_se = el(B,2,`j')
	}
	
	// Get HV deviance HV in validation and test data
	local cols : colfullnames E
	foreach m in newitems sameitems test {
		matrix E_`m' = E
		matrix colnames E_`m' = `=subinstr("`cols'", "_train", "_`m'", .)'
		`quietly' melogit y_`m' `varlist' || person:, asis from(E_`m') iterate(0)
		return scalar dev_`m' = -2 * e(ll)
	}
	
	// Get insample deviance from test data
	`quietly' melogit y_test `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_intest = -2 * e(ll)
	
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
	
	// Training data estimation
	`quietly' gsem (b_train <- c.se_train#c.M@1 `varlist'), variance(M@1)
	return scalar dev_train = -2 * e(ll)
	matrix B =  r(table)
	matrix E = e(b)
	if (e(converged) == 0) local ++converge_fails
	predict delta_hat_training, eta conditional(fixedonly)
	generate tausq = el(E, 1, colsof(E))
	
	// Get RMSE for HV
	rmse, predict_var(delta_hat_training)
	return scalar rmse_full = r(rmse_full)
	return scalar rmse_fix = r(rmse_fix)
	
	// Information criteria
	`quietly' estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)
	
	// Return parameter estimates
	local cols : colnames B
	local cols `=regexr("`cols'", "c.se_train#c.M", "drop")'
	local cols `=subinstr("`cols'", "x", "beta", .)'
	local cols `=subinstr("`cols'", "c.", "", .)'
	//local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
	//local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
	local cols `=subinstr("`cols'", "_cons", "beta1", 1)'
	local cols `=subinstr("`cols'", "_cons", "drop", 1)'
	local cols `=subinstr("`cols'", "_cons", "tausq", 1)'
	forvalues j = 1/`=colsof(B)' {
		local cname = word("`cols'", `j')
		if ("`cname'" != "drop") {
			return scalar `cname'_est = el(B, 1, `j')
			return scalar `cname'_se = el(B, 2, `j')
		}
	}
	
	// Holdout validation meta-regression
	foreach n in newitems sameitems test {
		`quietly' generate ll_`n' = ///
			log(normalden(b_`n', delta_hat, sqrt(tausq + se_`n'^2)))
		`quietly' summarize ll_`n'
		return scalar dev_`n' = -2*r(sum)
	}

	// In-sample meta-regression (test data only)
	`quietly' gsem (b_test <- c.se_test#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_intest = -2*e(ll)	
	
	// LOO CV meta-regression
	quietly generate ll_loo = .
	forvalues i = 1/`I' {
		`quietly' gsem (b_train <- c.se_train#c.M@1 `varlist') ///
			if item != `i', variance(M@1)
		if (e(converged) == 0) local ++converge_fails
		predict delta_hat_temp, eta conditional(fixedonly)
		matrix E_loo = e(b)
		local tau2 = el(E_loo, 1, colsof(E_loo))
		`quietly' replace ll_loo = log(normalden(delta_hat_temp, b_train, ///
			sqrt(`tau2' + se_train^2))) if item == `i'
		drop delta_hat_temp
	}
	`quietly' summarize ll_loo
	return scalar dev_loo = -2*r(sum)
	
	return scalar converge_fails = `converge_fails'
	
restore
end


// Simulation ------------------------------------------------------------------

// Import random number seeds
import delimited "random.txt", delimiter(tab) varnames(nonames) clear
keep in 1/510

// Set up simulation conditions, store in a matrix
egen condition = seq(), from(1) to(5)
generate double Vsq = 1.5^2
generate double Rsq = .6
replace Rsq = .3 if condition == 1
replace Rsq = .9 if condition == 3
generate double upsilon = sqrt(Rsq*Vsq)
generate double tau = sqrt(Vsq - upsilon^2)
generate double b = sqrt(upsilon^2 / 3)
generate nitems = 32
replace nitems = 16 if condition == 4
replace nitems = 64 if condition == 5
generate double sigma = 1.5
generate npersons = 500
generate seed = v`s'
mkmat *, matrix(SIM)

// Save a file describing each condition
keep in 1/5 
drop v* seed
save conditions, replace
clear

// Define models (predictors)
local m1 "x2-x3"
local m2 "x2-x4"
local m3 "x2-x5"

// Prepare to store results
postfile memhold_sim double(     ///
		seed                     ///
		condition                ///
		converge_fails           ///
		dev_meta                 ///
		dev_rasch                ///
	) using results_sim_`s', replace
postfile memhold_lltm double(    ///
		seed                     ///
		condition                ///
		model                    ///
		converge_fails           ///
		dev_intest     ///
		dev_test       ///
		dev_sameitems  ///
		dev_newitems   ///
		sigmasq_se     ///
		sigmasq_est    ///
		beta1_se       ///
		beta1_est      ///
		beta5_se       ///
		beta5_est      ///
		beta4_se       ///
		beta4_est      ///
		beta3_se       ///
		beta3_est      ///
		beta2_se       ///
		beta2_est      ///
		bic            ///
		aic            ///
		dev_train      ///
		rmse_full      ///
		rmse_fix       ///
	) using results_lltm_`s', replace
postfile memhold_twostage double(    ///
		seed                     ///
		condition                ///
		model                    ///
		converge_fails           ///
		dev_loo       ///
		dev_intest    ///
		dev_test      ///
		dev_sameitems ///
		dev_newitems  ///
		tausq_se      ///
		tausq_est     ///
		beta1_se      ///
		beta1_est     ///
		beta5_se      ///
		beta5_est     ///
		beta4_se      ///
		beta4_est     ///
		beta3_se      ///
		beta3_est     ///
		beta2_se      ///
		beta2_est     ///
		bic           ///
		aic           ///
		dev_train     ///
		rmse_full     ///
		rmse_fix      ///
	) using results_twostage_`s', replace

timer clear
timer on 4
forvalues i = 1/`=rowsof(SIM)' {

	display "$S_TIME: Starting `i' of `=rowsof(SIM)'"

	timer on 1
	sim_dataset, ///
		nitems(`=el(SIM, `i', colnumb(SIM, "nitems"))') ///
		npersons(`=el(SIM, `i', colnumb(SIM, "npersons"))') ///
		sigma(`=el(SIM, `i', colnumb(SIM, "sigma"))') ///
		tau(`=el(SIM, `i', colnumb(SIM, "tau"))') ///
		b(`=el(SIM, `i', colnumb(SIM, "b"))') ///
		seed(`=el(SIM, `i', colnumb(SIM, "seed"))')
	post memhold_sim ///
		(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
		(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
		(r(converge_fails))  ///
		(r(dev_meta)      )  ///
		(r(dev_rasch)     )
	timer off 1
	
	forvalues m = 1/3 {
	
		timer on 2
		capture noisily lltm_cv `m`m''
		post memhold_lltm ///
			(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
			(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
			(`m') ///
			(r(converge_fails)) ///
			(r(dev_intest   ))  ///
			(r(dev_test     ))  ///
			(r(dev_sameitems))  ///
			(r(dev_newitems ))  ///
			(r(sigmasq_se   ))  ///
			(r(sigmasq_est  ))  ///
			(r(beta1_se     ))  ///
			(r(beta1_est    ))  ///
			(r(beta5_se     ))  ///
			(r(beta5_est    ))  ///
			(r(beta4_se     ))  ///
			(r(beta4_est    ))  ///
			(r(beta3_se     ))  ///
			(r(beta3_est    ))  ///
			(r(beta2_se     ))  ///
			(r(beta2_est    ))  ///
			(r(bic          ))  ///
			(r(aic          ))  ///
			(r(dev_train    ))  ///
			(r(rmse_full)    )  ///
			(r(rmse_fix)     )
		timer off 2

		timer on 3
		capture noisily twostage_cv `m`m''
		post memhold_twostage ///
			(`=el(SIM, `i', colnumb(SIM, "seed"))') ///
			(`=el(SIM, `i', colnumb(SIM, "condition"))') ///
			(`m') ///
			(r(converge_fails)) ///
			(r(dev_loo      ))  ///
			(r(dev_intest   ))  ///
			(r(dev_test     ))  ///
			(r(dev_sameitems))  ///
			(r(dev_newitems ))  ///
			(r(tausq_se     ))  ///
			(r(tausq_est    ))  ///
			(r(beta1_se     ))  ///
			(r(beta1_est    ))  ///
			(r(beta5_se     ))  ///
			(r(beta5_est    ))  ///
			(r(beta4_se     ))  ///
			(r(beta4_est    ))  ///
			(r(beta3_se     ))  ///
			(r(beta3_est    ))  ///
			(r(beta2_se     ))  ///
			(r(beta2_est    ))  ///
			(r(bic          ))  ///
			(r(aic          ))  ///
			(r(dev_train    ))  ///
			(r(rmse_full)    )  ///
			(r(rmse_fix)     )
		timer off 3
		
	}

}

timer off 4
timer list

postclose memhold_sim
postclose memhold_twostage
postclose memhold_lltm

// Save an example dataset
// clear
// sim_dataset
// save example, replace
