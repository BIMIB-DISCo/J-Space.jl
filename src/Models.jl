#using DataFrames


# π are basefrequencies and they do not need to  sum to 1
#example
#π=[3.,4.,5.,6.]

# α is the transitional substitutions rate or
# the transitional over transversal rate ratio

# β is the transversional substitution rate

# αʳ is the purines (A<->G) transitional substitutions rate
# αʸ is the pyrimidines(Cz->T) transitional substitutions rate

#βᶜᵍ⁻ᵃᵗ   transversal switching with the paired base C<->G A<->T
#βᵃᶜ⁻ᵍᵗ   transversal rate A<->C G<->T

#note: the  previous two could be due to error corrections
#(biological origin not specified in the article)

#βᵗʳ thymine from/to purine
#βᶜʳ cytosine from/to purine

#example
#par=IdDict( "alpha" => 0.5, "beta" => 0.3)





#makes equal stationary frequencies
function setequalfreqs()
    return(repeat([0.25],1,4))
end

#generate normalized rate matrix Q for different models
#with main diagonal set to 0 instead of -1
#π are the stationary base frequencies ordered as [:A, :C, :G, :T]
#par is a dictionary with rates (each model in its own way)
function Q(model::String, par::IdDict)

    #println("============")
    #println("model $model")

    preferred_order = [:A, :C, :G, :T]
    ord=indexin(preferred_order,[:T, :C, :A, :G])

    a=b=c=d=e=f=a1=b1=c1=d1=e1=f1=-1;


	#######
	#check input parameters

	#check π
	if ( model == "F81" || model == "HKY" || model == "TrN" || model == "K81uf")
		if "pi" ∉ keys(par)
			println("Warning -> Stationary frequencies not provided will all be set to 1/4 ")
			π=setequalfreqs()
		else
			π=par["pi"]
		end
    else
        π=setequalfreqs()
	end
	#println("par: ",par)
	#check α
	if ( model == "K80" || model == "HKY" || model == "TrN93ef" || model == "TrN" )
		if "alpha" ∉ keys(par)
			println("Error -> Parameter α for the model not provided")
			return "Error"
		else
			α=par["alpha"]
		end
	elseif ( model == "JC69" || model == "F81" || model == "K81" || model == "K81uf" )
		if "alpha" ∉ keys(par)
			println("Warning -> Parameter α for the model not provided wil be set to 1")
			α=1.0
		else
			α=par["alpha"]
		end
	end

	#check α2
	if (model == "TrN93ef" || model == "TrN" )
		if "alpha2" ∉ keys(par)
			println("Error -> Parameter α2 for the model not provided")
			return "Error"
		else
			α2=par["alpha2"]
		end
	end

	#check β
	if ( model == "K81" || model == "K81uf" )
		if "beta" ∉ keys(par)
			println("Error -> Parameter β for the model not provided")
			return "Error"
		else
			β=par["beta"]
		end
	elseif ( model == "K80" || model == "HKY" || model == "TrN93ef" || model == "TrN" )
		if "beta" ∉ keys(par)
			println("Warning -> Parameter β for the model set to 1")
			β=1.0
		else
			β=par["beta"]
		end
	end

	#check β2
	if (model == "K81" || model == "K81uf")
		if "beta2" ∉ keys(par)
			println("Error -> Parameter β2 for the model not provided")
			return "Error"
		else
			β2=par["beta2"]
        end
	end
	# end checking
	#######

	#order π following a fixed pattern
	π = π[ord]

    if ( model == "JC69" || model == "F81" )
        if ( model == "JC69" )
            π=setequalfreqs()
        end

        a=b=c=d=e=f=a1=b1=c1=d1=e1=f1=α;

    elseif ( model == "K80" || model == "HKY" )
        if ( model == "K80" )
            π=setequalfreqs()
		end

		a=f=α;
		b=c=d=e=β
        #("beta" ∈ keys(par)) ? b=c=d=e=par["beta"] : b=c=d=e=1;  #check and delete

    elseif ( model == "TrN93ef" || model == "TrN" )
        #Tn93: Tamura, K. and Nei, M. (1993) Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Mol. Bio. Evol. 10:512-526.
		if ( model == "TrN93ef" )
            π=setequalfreqs()
        end

		αʸ=α;
		αʳ=α2;

		a=αʸ;
		b=c=d=e=β;
        #("beta" ∈ keys(par)) ? b=c=d=e=par["beta"] : b=c=d=e=1;  #check and delete
        f=αʳ;

	elseif ( model == "K81" || model == "K81uf" )
		#K81: Kimura, M. (1981) Estimation of evolutionary distances between homologous nucleotide sequences. Proc. Natl. Acad. Sci. USA. 78:454-458.
        if ( model == "K81" )
            π=setequalfreqs()
        end

		βᶜᵍ⁻ᵃᵗ = par["beta"]
		βᵃᶜ⁻ᵍᵗ = par["beta2"]

		a=f=α
		#("alpha" ∈ keys(par)) ? a=f=par["alpha"] : a=f=1;  #check and delete
        b=e=βᶜᵍ⁻ᵃᵗ;
        c=d=βᵃᶜ⁻ᵍᵗ;
    else
        println("model ", model, " not found")
    end

    if ( model != "ERR" || model != "WRONG")

			#make matrix symmetrical
			a1=a; b1=b; c1=c; d1=d; e1=e; f1=f;

            #println("test for errors")

            err_string="[MODEL] $model INT:makeQmatrix; " *
                "DNA substitution parameter routine error -"*
                "one or more values = -1" #alighieri said

            @assert !( a==-1 || b==-1 || c==-1 || d==-1 ||
                e==-1 || f==-1 || a1==-1 || b1==-1 ||
                c1==-1 || d1==-1 || e1==-1 || f1==-1 ) err_string


			#println("no errors")
	end

    #set entries of Q matrix.
    Q=DataFrame(
        [0  a  b  c ;
        a1 0  d  e ;
        b1 d1 0  f ;
        c1 e1 f1 0 ],
        [:T, :C, :A, :G])

    sum=0.
    if ( model != "ERR" || model != "WRONG" )
		#rescale stationary frequencies, to make certain they sum to 1.0
		for i ∈ 1:4
            sum += π[i]
        end
		if ( sum ≠ 1 )
            for i ∈ 1:4
                π[i] /= sum
            end
        end

		#multiply entries by stationary frequencies
		for i ∈ 1:4, j ∈ 1:4
            Q[i,j] *= π[j]; #base freqs are same along a column
        end
    end

    #rescale, so branch lengths are in terms of expected number
    #of substitutions per site
    #calculate scale factor and note that Q_ii = 0
    #=scaler = 0.0
    for i ∈ 1:4, j ∈ 1:4
        scaler += π[i] * Q[i,j]
    end

    scaler = 1.0 / scaler
    for i ∈ 1:4, j ∈ 1:4
        Q[i,j] *= scaler
    end=#



    #set diagonal of matrix
    # for i ∈ 1:4
    #     sum=0;
    #     for j ∈ 1:4
    #         sum += Q[i,j]
    #     end
    #     Q[i,i]=-sum
    # end


    Q=Q[ord,ord]
    #println("Q = $Q")

    return Q;
end
#=
#Q("abc")
par=IdDict()
println(Q("JC69", par))
println(Q("K80", par))

par=IdDict( "alpha" => 0.5)
println(Q("JC69", par))
println(Q("K80", par))
println(Q("K81", par))

par=IdDict( "alpha" => 0.5, "beta" => 0.3)
println(Q("JC69", par))
println(Q("K80", par))
println(Q("K81", par))

par=IdDict( "beta" => 0.3, "beta2" => 0.3)
println(Q("JC69", par))
println(Q("K80", par))
println(Q("K81", par))

par=IdDict( "alpha" => 0.3, "alpha2" => 0.3)
println(Q("TrN93ef", par))
=#
