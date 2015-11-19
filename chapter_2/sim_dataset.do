capture program drop sim_dataset
program define sim_dataset
	
	syntax [, i_factor(integer 3) p_factor(integer 100) ///
		tau(real .3) sigma(real 1) ]
	
	local beta1 = -.5
	local beta2 = .5
	local beta3 = .5
	local beta4 = .5
	local beta5 = -.5
	local gamma1 = .5
	local gamma2 = .5
	
	clear
	
	// Make item 8-item set (all x* combinations)
	quietly: set obs 2
	generate x1 = 1
	forvalues i = 2/4 {
		egen x`i' = seq(), from(0) to(1)
	}
	fillin x*
	drop _fillin
	
	// Expand to full item set. Sim residuals and calc difficulty.
	quietly: expand `=`i_factor''
	quietly: expand 2, generate(new_item)
	generate item = _n
	generate epsilon = rnormal(0, `tau')
	generate xB = `beta1'*x1 + `beta2'*x2 + `beta3'*x3 + `beta4'*x4 + `beta5'*x2*x3
	generate delta = xB + epsilon
	
	if `p_factor' > 0 {
	
		// Make 4-person set (all w* combinations)
		quietly: expandcl 4, cluster(x1) generate(w_temp)
		quietly: recode w_temp (1 2 = 0) (else = 1), generate(w1)
		quietly: recode w_temp (1 3 = 0) (else = 1), generate(w2)
		quietly: expand 2, generate(new_person)
		quietly: expandcl `p_factor', cluster(new_person w_temp) generate(person)
		
		// Add person residuals, calculate fixed and full ability
		quietly: generate zeta_temp = rnormal(0, `sigma') if item == 1
		quietly: bysort person: egen zeta = mean(zeta_temp)
		generate wC = `gamma1'*w1 + `gamma2'*w2
		generate theta = wC + zeta
		
		drop *_temp
		
		// Simulate response
		generate yhat = invlogit(theta - delta)
		generate y = rbinomial(1, yhat)
		
	}	

end
