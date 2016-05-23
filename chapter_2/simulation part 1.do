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
	generate epsilon = rnormal(0, `tau')
	generate epsilon_newitems = rnormal(0, `tau')
	
	// Expand dataset to P persons and generate abilities
	quietly: expandcl `npersons', cluster(temp_one) generate(person)
	quietly: generate temp_theta = rnormal(0, `sigma') if item == 1
	quietly: generate temp_theta_new = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen theta = mean(temp_theta)
	quietly: bysort person: egen theta_new = mean(temp_theta_new)

	// Simulate responses
	generate y = rbinomial(1, invlogit(theta - xB - epsilon))
	generate y_newpersons = rbinomial(1, invlogit(theta_new - xB - epsilon))
	generate y_newitems = rbinomial(1, invlogit(theta - xB - epsilon_new))
	generate y_newboth = rbinomial(1, invlogit(theta_new - xB - epsilon_new))
	
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
	quietly generate `squared_error' = (-1*`delta_hat' - (xB + epsilon))^2
	quietly summarize `squared_error'
	return scalar rmse = sqrt(r(mean))
	
end


// Function to conduct LLTM analysis -------------------------------------------

capture program drop lltm_cv
program define lltm_cv, rclass
preserve

	syntax varlist (fv)

	local converge_fails = 0
	
	// Estimate model on training data with new persons
	quietly melogit y_newpersons `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	matrix B_newpersons = e(b)
	local cols : colfullnames B_newpersons 
	matrix colnames B_newpersons = `=subinstr("`cols'", "_newpersons", "", .)'
	rmse, predict_opts(xb)
	return scalar rmse_newpersons = r(rmse)
	
	// Estimate model on training data with new items
	quietly melogit y_newitems `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	matrix B_newitems = e(b)
	local cols : colfullnames B_newitems 
	matrix colnames B_newitems = `=subinstr("`cols'", "_newitems", "", .)'
	rmse, predict_opts(xb)
	return scalar rmse_newitems = r(rmse)
	
	// Estimate model on main dataset
	quietly melogit y `varlist' || person:
	if (e(converged) == 0) local ++converge_fails
	return scalar dev_insample = -2 * e(ll)
	rmse, predict_opts(xb)
	return scalar rmse_insample = r(rmse)
	
	// Get IC for insample fit
	quietly estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)
	
	// Apply new person hold out fit to main dataset
	quietly melogit y `varlist' || person:, asis from(B_newpersons) iterate(0)
	return scalar dev_newpersons = -2 * e(ll)
	
	// Apply new item hold out fit to main dataset
	quietly melogit y `varlist' || person:, asis from(B_newitems) iterate(0)
	return scalar dev_newitems = -2 * e(ll)
	
	return scalar converge_fails = `converge_fails'
	
restore
end


// Function to conduct two-stage analysis --------------------------------------

capture program drop twostage_cv
program define twostage_cv, rclass
preserve

	syntax varlist (fv)
	
	local converge_fails = 0
	quietly summarize item
	local I = r(max)

	// Get Rasch difficulties
	quietly melogit y ibn.item, noconstant || person:
	if (e(converged) == 0) local ++converge_fails
	matrix A = r(table)
	matrix B = A[1..2, 1..`I']'

	// Get Rasch difficulties from new item training data
	quietly melogit y_newitems ibn.item, noconstant || person:
	if (e(converged) == 0) local ++converge_fails
	matrix A = r(table)
	matrix B_newitems = A[1..2, 1..`I']'

	// Get Rasch difficulties from new person training data
	quietly melogit y_newpersons ibn.item, noconstant || person:
	if (e(converged) == 0) local ++converge_fails
	matrix A = r(table)
	matrix B_newpersons = A[1..2, 1..`I']'
	
	// Make dataset of Rasch difficultes
	quietly keep if person == 1
	keep item x2-x4 xB epsilon
	sort item
	svmat B_newitems, names(col)
	rename (b se) (b_newitems se_newitems)
	svmat B_newpersons, names(col)
	rename (b se) (b_newpersons se_newpersons)
	svmat B, names(col)
	
	// In-sample meta-regression
	quietly gsem (b <- c.se#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	local dev_insample = -2*e(ll)
	return scalar dev_insample = `dev_insample'
	rmse, predict_opts(eta conditional(fixedonly))
	return scalar rmse_insample = r(rmse)
	
	// TRMSEA
	generate ll_saturated = log(normalden(0, 0, se))
	quietly summarize ll_saturated
	local dev_saturated = -2*r(sum)
	local df = e(N) - e(rank)
	return scalar rmsea = ///
		sqrt(1/(e(N)*`df') * max(`dev_insample' - `dev_saturated' - `df', 0))
	
	// Get IC
	quietly estat ic
	return scalar aic = el(r(S), 1, 5)
	return scalar bic = el(r(S), 1, 6)

	// Holdout validation meta-regression with new items
	quietly gsem (b_newitems <- c.se_newitems#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	predict delta_hat, eta conditional(fixedonly)
	matrix E = e(b)
	local tau2 = el(E, 1, colsof(E))
	generate ll_newitems = ///
		log(normalden(b, delta_hat, sqrt(`tau2' + se_newitems^2)))
	quietly summarize ll_newitems
	return scalar dev_newitems = -2*r(sum)
	rmse, predict_var(delta_hat)
	return scalar rmse_newitems = r(rmse)
	drop delta_hat

	// Holdout validation meta-regression with new items
	quietly gsem (b_newpersons <- c.se_newpersons#c.M@1 `varlist'), variance(M@1)
	if (e(converged) == 0) local ++converge_fails
	predict delta_hat, eta conditional(fixedonly)
	matrix E = e(b)
	local tau2 = el(E, 1, colsof(E))
	generate ll_newpersons = ///
		log(normalden(b, delta_hat, sqrt(`tau2' + se_newpersons^2)))
	quietly summarize ll_newpersons
	return scalar dev_newpersons = -2*r(sum)
	rmse, predict_var(delta_hat)
	return scalar rmse_newpersons = r(rmse)
	drop delta_hat
	
	// LOO CV meta-regression
	quietly generate ll_loo = .
	quietly generate delta_hat = .
	forvalues i = 1/`I' {
		quietly gsem (b <- c.se#c.M@1 `varlist') if item != `i', variance(M@1)
		if (e(converged) == 0) local ++converge_fails
		predict delta_hat_temp, eta conditional(fixedonly)
		matrix E = e(b)
		local tau2 = el(E, 1, colsof(E))
		quietly replace ll_loo = ///
			log(normalden(b, delta_hat_temp, sqrt(`tau2' + se^2))) if item == `i'
		quietly replace delta_hat = delta_hat_temp if item == `i'
		drop delta_hat_temp
	}
	quietly summarize ll_loo
	return scalar dev_loo = -2*r(sum)
	rmse, predict_var(delta_hat)
	return scalar rmse_loo = r(rmse)
	
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
postfile memhold_lltm double(seed) condition double(tau) nitems npersons model /// 
	double(dev_insample aic bic dev_newitems dev_newpersons ///
	rmse_insample rmse_newitems rmse_newpersons) ///
	using results_lltm_`s', replace
postfile memhold_twostage double(seed) condition double(tau) nitems npersons ///
	model ///
	double(dev_insample aic bic dev_newitems dev_newpersons dev_loo ///
	rmse_insample rmse_newitems rmse_newpersons rmse_loo rmsea) ///
	using results_twostage_`s', replace
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
				(r(dev_insample)) ///
				(r(aic)) ///
				(r(bic)) ///
				(r(dev_newitems)) ///
				(r(dev_newpersons)) ///
				(r(rmse_insample)) ///
				(r(rmse_newitems)) ///
				(r(rmse_newpersons))
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
				(r(dev_insample)) ///
				(r(aic)) ///
				(r(bic)) ///
				(r(dev_newitems)) ///
				(r(dev_newpersons)) ///
				(r(dev_loo)) ///
				(r(rmse_insample)) ///
				(r(rmse_newitems)) ///
				(r(rmse_newpersons)) ///
				(r(rmse_loo)) ///
				(r(rmsea))
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
