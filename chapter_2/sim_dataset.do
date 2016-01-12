capture program drop sim_dataset
program define sim_dataset
	
	syntax [, nitems(integer 32) npersons(integer 500) tau(real .5) sigma(real 1) ]
	
	local beta1 = -.375
	local beta2 = 1
	local beta3 = 1
	local beta4 = -1
	local beta5 = -.5
	local gamma1 = -.5
	local gamma2 = .5
	
	clear
	
	// Make base dataset (each unique combination of person and item fixed parts)
	quietly: set obs 2
	forvalues i = 1/2 {
		egen w`i' = seq(), from(0) to(1)
	}
	generate x1 = 1
	forvalues i = 2/4 {
		egen x`i' = seq(), from(0) to(1)
	}
	fillin w* x*

	// Fixed part of person and items
	generate wC = `gamma1'*w1 + `gamma2'*w2
	generate xB = `beta1'*x1 + `beta2'*x2 + `beta3'*x3 + `beta4'*x4 + `beta5'*x2*x3
	
	// Expand to full dataset
	egen temp_item_base = group(x1-x4)
	egen temp_person_base = group(w1-w2)
	quietly: expandcl `=`nitems'/8', cluster(temp_item_base) generate(item)
	quietly: expandcl `=`npersons'/4', cluster(temp_person_base) generate(person)
	
	// Simulate item residuals for in and out of sample
	quietly: generate temp_epsilon_in = rnormal(0, `tau') if person == 1
	quietly: generate temp_epsilon_out = rnormal(0, `tau') if person == 1
	quietly: bysort item: egen epsilon_in = mean(temp_epsilon_in)
	quietly: bysort item: egen epsilon_out = mean(temp_epsilon_out)
	
	// Simulate person residuals for in and out of sample
	quietly: generate temp_zeta_in = rnormal(0, `sigma') if item == 1
	quietly: generate temp_zeta_out = rnormal(0, `sigma') if item == 1
	quietly: bysort person: egen zeta_in = mean(temp_zeta_in)
	quietly: bysort person: egen zeta_out = mean(temp_zeta_out)
	
	// Simulate responses
	generate y = rbinomial(1, invlogit(wC + zeta_in - xB - epsilon_in))
	generate y_newpersons = rbinomial(1, invlogit(wC + zeta_out - xB - epsilon_in))
	generate y_newitems = rbinomial(1, invlogit(wC + zeta_in - xB - epsilon_out))
	generate y_newboth = rbinomial(1, invlogit(wC + zeta_out - xB - epsilon_out))
	
	drop _fillin temp_*

end

