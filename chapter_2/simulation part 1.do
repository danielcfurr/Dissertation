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
	local beta5 = `b'
	
	// Use seed if one is provided
	clear
	if(`seed' > -1) set seed `seed'
	
	// Set up the items for all subsets but same items validation
	quietly set obs 3
	quietly generate subset = "training"
	quietly replace subset = "newitems" in 2
	quietly replace subset = "evaluation" in 3
	quietly expandcl `nitems', cluster(subset) generate(item)
	generate x1 = 1
	forvalues i = 2/4 {
		generate x`i' = rnormal()
	}
	generate xB = ///
		`beta1'*x1 + `beta2'*x2 + `beta3'*x3 + `beta4'*x4 + `beta5'*x2*x3
	generate epsilon = rnormal(0, `tau')

	// Add in items from same items validation by copying from training subset
	quietly expand 2 if subset == "training", generate(temp_item)
	quietly replace subset = "sameitems" if temp_item
	
	// Expand dataset to P persons and generate abilities
	quietly generate temp_one = 1
	quietly expandcl `npersons', cluster(temp_one) generate(temp_person)
	quietly generate pick = temp_person == 1
	egen person = group(subset temp_person)
	egen temp_tag = tag(subset person)
	quietly generate temp_theta = rnormal(0, `sigma') if temp_tag
	quietly bysort person: egen theta = mean(temp_theta)
	sort person item
	
	// Simulate responses
	generate y = rbinomial(1, invlogit(theta - xB - epsilon))

	// Set up meta-analysis dataset
	local converge_fails = 0
	quietly generate delta = .
	quietly generate delta_se = .
	foreach s in training newitems sameitems evaluation {
		quietly melogit y ibn.item if subset == "`s'", noconstant || person:
		if (e(converged) == 0) local ++converge_fails
		// Return training sample Rasch deviance
		if ("`s'" == "training") return scalar dev_rasch = -2*e(ll)
		quietly predict temp_delta if subset == "`s'", xb
		quietly replace delta = temp_delta if subset == "`s'"
		quietly predict temp_delta_se if subset == "`s'", stdp
		quietly replace delta_se = temp_delta_se if subset == "`s'"
		drop temp_*
	}
	
	// Return training sample "Rasch" meta-regression deviance
	quietly generate temp_ll = log(normalden(0, 0, delta_se))
	quietly summarize temp_ll if subset == "training"
	return scalar dev_meta = -2 * r(sum)
			
	return scalar converge_fails = `converge_fails'
	
	drop temp_*
	
end


// Function to calculate root mean squared error for item predictions ----------

capture program drop rmse
program define rmse, rclass

	syntax ,  [predict_opts(string) predict_var(varname)]
	
	tempvar delta_hat squared_error
	
	if (e(cmd) == "meglm") local if "if pick"
	
	// note: delta_hat is easiness, so sign is reversed
	if ("`predict_var'" == "") {
		quietly predict `delta_hat' `if', `predict_opts'
	}
	else {
		quietly generate `delta_hat' = `predict_var' `if'
	}
	
	quietly generate `squared_error' = (-1*`delta_hat' - xB + epsilon)^2
	quietly summarize `squared_error' if subset == "evaluation"
	return scalar rmse_full = sqrt(r(mean))
	
	quietly replace `squared_error' = (-1*`delta_hat' - xB)^2
	quietly summarize `squared_error' if subset == "evaluation"
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
	`quietly' melogit y `varlist' if subset == "training" || person:
	return scalar dev_training = -2 * e(ll)
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
	local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
	local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
	local cols `=subinstr("`cols'", "_cons", "beta1", 1)'
	local cols `=subinstr("`cols'", "_cons", "sigmasq", 1)'
	forvalues j = 1/`=colsof(B)' {
		local cname = word("`cols'", `j')
		return scalar `cname'_est = el(B,1,`j')
		return scalar `cname'_se = el(B,2,`j')
	}
	
	// Get HV deviance HV in validation and evaluation data
	local cols : colfullnames E
	foreach s in newitems sameitems evaluation {
		matrix E_`s' = E
		matrix colnames E_`s' = `=subinstr("`cols'", "_training", "_`s'", .)'
		`quietly' melogit y `varlist' if subset == "`s'" || person:, ///
			asis from(E_`s') iterate(0)
		return scalar dev_`s' = -2 * e(ll)
	}
	
	// Get insample deviance from evaluation data
	`quietly' melogit y `varlist' if subset == "evaluation" || person:
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_ineval = -2 * e(ll)
	
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
	`quietly' keep if pick
	local I = _N / 4
	
	// Training data estimation
	`quietly' gsem (delta <- c.delta_se#c.M@1 `varlist') ///
		if subset == "training", variance(M@1)
	return scalar dev_training = -2 * e(ll)
	matrix B =  r(table)
	matrix E = e(b)
	if (e(converged) == 0) local ++converge_fails
	quietly predict delta_hat, eta conditional(fixedonly)
	quietly generate tausq = el(E, 1, colsof(E))
	
	// Get RMSE for HV
	rmse, predict_var(delta_hat)
	return scalar rmse_full = r(rmse_full)
	return scalar rmse_fix = r(rmse_fix)
	
	// Information criteria
	`quietly' estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)
	
	// Return parameter estimates
	local cols : colnames B
	local cols `=regexr("`cols'", "c.delta_se#c.M", "drop")'
	local cols `=subinstr("`cols'", "x", "beta", .)'
	local cols `=subinstr("`cols'", "c.", "", .)'
	local cols `=subinstr("`cols'", "beta2#beta3", "beta5", .)'
	local cols `=subinstr("`cols'", "beta2#beta4", "beta6", .)'
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
	`quietly' generate ll = ///
		log(normalden(delta, delta_hat, sqrt(tausq + delta_se^2)))
	foreach s in newitems sameitems evaluation {
		`quietly' summarize ll if subset == "`s'"
		return scalar dev_`s' = -2*r(sum)
	}

	// In-sample meta-regression (evaluation data only)
	`quietly' gsem (delta <- c.delta_se#c.M@1 `varlist') ///
		if subset == "evaluation", variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_ineval = -2*e(ll)	
	
	// LOO CV meta-regression
	quietly generate ll_loo = .
	quietly generate delta_hat_loo = .
	`quietly' levels item if subset == "training", local(i_list)
	foreach i of local i_list {
		`quietly' gsem (delta <- c.delta_se#c.M@1 `varlist') ///
			if item != `i' & subset == "training", variance(M@1)
		if (e(converged) == 0) local ++converge_fails
		quietly predict temp_delta_hat if item == `i' & ///
			subset == "training", eta conditional(fixedonly)
		quietly replace delta_hat_loo = temp_delta_hat ///
			if item == `i' & subset == "training"
		matrix E_loo = e(b)
		local tau2 = el(E_loo, 1, colsof(E_loo))
		quietly replace ll_loo = log(normalden(delta_hat_loo, delta, ///
			sqrt(`tau2' + delta_se^2))) if item == `i' & subset == "training"
		drop temp_*
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
generate double b = sqrt(upsilon^2 / 4)
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
local m1 "x2-x4"
local m2 "x2-x4 c.x2#c.x3"
local m3 "x2-x4 c.x2#c.x3 c.x2#c.x4"

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
		dev_ineval     ///
		dev_evaluation       ///
		dev_sameitems  ///
		dev_newitems   ///
		sigmasq_se     ///
		sigmasq_est    ///
		beta1_se       ///
		beta1_est      ///
		beta6_se       ///
		beta6_est      ///
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
		dev_training      ///
		rmse_full      ///
		rmse_fix       ///
	) using results_lltm_`s', replace
postfile memhold_twostage double(    ///
		seed                     ///
		condition                ///
		model                    ///
		converge_fails           ///
		dev_loo       ///
		dev_ineval    ///
		dev_evaluation      ///
		dev_sameitems ///
		dev_newitems  ///
		tausq_se      ///
		tausq_est     ///
		beta1_se      ///
		beta1_est     ///
		beta6_se       ///
		beta6_est      ///
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
		dev_training     ///
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
			(r(dev_ineval   ))  ///
			(r(dev_evaluation     ))  ///
			(r(dev_sameitems))  ///
			(r(dev_newitems ))  ///
			(r(sigmasq_se   ))  ///
			(r(sigmasq_est  ))  ///
			(r(beta1_se     ))  ///
			(r(beta1_est    ))  ///
			(r(beta6_se     ))  ///
			(r(beta6_est    ))  ///
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
			(r(dev_training    ))  ///
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
			(r(dev_ineval   ))  ///
			(r(dev_evaluation     ))  ///
			(r(dev_sameitems))  ///
			(r(dev_newitems ))  ///
			(r(tausq_se     ))  ///
			(r(tausq_est    ))  ///
			(r(beta1_se     ))  ///
			(r(beta1_est    ))  ///
			(r(beta6_se     ))  ///
			(r(beta6_est    ))  ///
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
			(r(dev_training    ))  ///
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
