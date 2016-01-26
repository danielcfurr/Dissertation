// Run in pseudo-parallel by opening multiple instances of Stata and changing
// 'instance' for each.
local instance = 1
local replications = 20

run sim_dataset.do

capture program drop cv_over_persons
program define cv_over_persons, rclass

	syntax varlist(fv) [, pcheck(varname fv max=1) *]
	
	melogit `varlist' if !new_person & !new_item || person:, `options'
	matrix define B = e(b)
	return scalar p_converged = e(converged) 
	return scalar p_dev_00 = -2*e(ll)
	
	matrix define T = r(table)
	return matrix table = T, copy
	if "`pcheck'" == "" {
		return scalar p_pcheck = .
	}
	else {
		local Tcol = "y:" + "`pcheck'"
		return scalar p_pcheck = T[rownumb(T, "pvalue"), colnumb(T, "`Tcol'")]
	}
	
	estat ic
	return scalar p_aic = el(r(S), 1, 5)
	
	melogit `varlist' if !new_item & new_person || person:, `options' ///
		asis from(B) iterate(0)
	return scalar p_dev_01 = -2*e(ll)

	melogit `varlist' if new_item & !new_person || person:, `options' ///
		asis from(B) iterate(0)
	return scalar p_dev_10 = -2*e(ll)

	melogit `varlist' if new_item & new_person || person:, `options' ///
		asis from(B) iterate(0)
	return scalar p_dev_11 = -2*e(ll)
	
end


capture program drop cv_over_items
program define cv_over_items, rclass

	syntax varlist(fv) [, *]
	tempvar rawscore one rawscore_in_train tag_item tag_person usable
	
	quietly bysort person new_person new_item: egen `rawscore' = total(y)
	quietly generate `one' = 1 if !new_person & !new_item
	quietly egen `rawscore_in_train' = mean(`one'), by(`rawscore')
	quietly replace `rawscore_in_train' = 0 if `rawscore_in_train' == .

	quietly egen `tag_item' = tag(item new_person)
	quietly summarize item if !new_person & !new_item & `tag_item'
	global maxscore = r(N)

	quietly egen `tag_person' = tag(person new_item)
	quietly summarize person if !new_person & !new_item & `tag_person'
	global P = r(N)

	quietly generate `usable' = `rawscore_in_train' & !inlist(`rawscore', 0, $maxscore)
	quietly summarize person if `usable' & `tag_person' & !new_person & !new_item
	return scalar proploss_00 = 1 - r(N) / $P
	quietly summarize person if `usable' & `tag_person' & new_person & new_item
	return scalar proploss_11 = 1 - r(N) / $P

	melogit `varlist' i.`rawscore' if `usable' & !new_person & !new_item || item:, `options'
	return scalar i_converged = e(converged) 
	return scalar i_dev_00 = -2*e(ll)
	matrix B = e(b)
	estat ic
	return scalar i_aic = el(r(S), 1, 5)
	
	melogit `varlist' i.`rawscore' if `usable' & new_person & new_item || item: , ///
		asis from(B, skip) iterate(0) `options'
	return scalar i_dev_11 = -2*e(ll)

end



clear
set obs 4
generate double tau = 0
replace tau = .1 in 2
replace tau = .3 in 3
replace tau = .5 in 4
generate pfactor = 250
replace pfactor = 75 in 1
generate ifactor = 16
replace ifactor = 4 in 1
fillin tau pfactor ifactor
drop _fillin

generate condition = _n
//generate state = ""
expand `replications'

local N = _N
local M = 3

local p_model_1 "y w1 w2 x2-x4"
local p_model_2 "y w1 w2 x2-x4 c.x2#c.x3, pcheck(c.x2#c.x3)"
local p_model_3 "y w1 w2 x2-x4 c.x2#c.x3 c.x3#c.x4, pcheck(c.x3#c.x4)"
local p_vars "p_dev_00 p_dev_10 p_dev_01 p_dev_11 p_pcheck p_aic p_converged"

local i_model_1 "y x2-x4"
local i_model_2 "y w1 w2 x2-x4 c.x2#c.x3"
local i_model_3 "y w1 w2 x2-x4 c.x2#c.x3 c.x3#c.x4"
local i_vars "i_dev_00 i_dev_11 proploss_00 proploss_11 i_aic i_converged"

forvalues m = 1/`M' {
	foreach v of local p_vars {
		generate double `v'_`m' = .
	}
	foreach v of local i_vars {
		generate double `v'_`m' = .
	}
}

timer clear
timer on 1

forvalues n = 1/`N' {
	local tau = tau[`n']
	local ifactor = ifactor[`n']
	local pfactor = pfactor[`n']
	display "$S_TIME: `n' of `N'"
	display "  tau = `tau', ifactor = `ifactor', pfactor = `pfactor'"
	preserve
		//local state = c(rngstate_mt64)
		sim_dataset, tau(`tau') p_factor(`pfactor') i_factor(`ifactor')
		forvalues m = 1/`M' {
			quietly cv_over_persons `p_model_`m''
			foreach v of local p_vars {
				local `v'_`m' = r(`v')
			}
			if `tau' > 0 {
				quietly cv_over_items `i_model_`m''
				foreach v of local i_vars {
					local `v'_`m' = r(`v')
				}
			}
		}
	restore
	forvalues m = 1/3 {
		foreach v of local p_vars {
			quietly replace `v'_`m' = ``v'_`m'' in `n'
		}
		if `tau' > 0 {
			foreach v of local i_vars {
				quietly replace `v'_`m' = ``v'_`m'' in `n'
			}
		}
	}
	//quietly: replace state = "`state'"
}

timer off 1
timer list

save "sim_`instance'.dta", replace
