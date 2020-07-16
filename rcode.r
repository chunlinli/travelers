library(ROCR)
library(stringi)

#################################################
# preprocessing data
# remove missing values, remove ni.age >= 100

data.train = read.csv("pp_train.csv")
data.train = na.omit(data.train)
data.train = data.train[data.train$ni.age < 100,]
data.train = data.train[,1:17]

#################################################
# this function calculates err

misclassification = function(obs,pred){
	temp = ifelse( obs == pred, 0, ifelse(obs == 1, 5, 1) )
	mean(temp)
}

#################################################
# step 1: L0 penalization selection of main effects
# this function returns optimal penalization parameter 2:log(n) by CV

L0Selection.tune = function(cv.k=100, cv.ratio=0.2, 
					lambda.step=0.1, lambda.l=1, lambda.u=7.2){
	n.data = nrow(data.train)

	lambda = seq(from = lambda.l, to = lambda.u, by = lambda.step)
	cv.err = numeric(length(lambda))
	cv.auc = numeric(length(lambda))

	for (j in 1:length(lambda)){
	
		err = numeric(cv.k)
		cstat = numeric(cv.k)
		
		for (i in 1:cv.k) {
			train = sample(1:n.data)[1:round(n.data*cv.ratio)]
			logit = glm(cancel ~ factor(claim.ind) + sales.channel + credit + 
			factor(zip.code) + n.adults + n.children + tenure + 
			len.at.res + ni.age + ni.marital.status + premium + 
			ni.gender + house.color + year + dwelling.type + coverage.type, 
			data = data.train[train,], family = binomial(link = "logit") )

			temp = step(logit, direction = "both", trace = 0, k = lambda[j])
			logit = eval(temp$call)
			pred = predict(logit, data.train[-train,], type = "response")

			input = prediction(pred, data.train[-train,]$cancel)
			auc = performance(input, "auc")
			cstat[i] = unlist(auc@y.values)

			pred = ifelse(pred > 1/6,1,0)
			err[i] = with(data.train[-train,], misclassification(cancel,pred))
		}

		cv.err[j] = mean(err)
		cv.auc[j] = mean(cstat)
	}

	idx = which.min(cv.err - cv.auc)
	return(lambda[idx])
}

# this function returns the counts of appearances of each predictors

L0Selection.count = function(cv.k = 500, cv.ratio = 0.2, lambda = 2){

	c = numeric(16)
	names(c) = c("factor(claim.ind)", "sales.channel", "credit", "factor(zip.code)", 
			"n.adults", "n.children", "tenure",  "len.at.res", 
			"ni.age","ni.marital.status", "premium", "ni.gender", "house.color", 
			"year", "dwelling.type", "coverage.type")

	for (i in 1:cv.k) {
	
		train = sample(1:n.data)[1:round(n.data*cv.ratio)]
		logit = glm(cancel ~ factor(claim.ind) + sales.channel + credit + 
			factor(zip.code) + n.adults + n.children + tenure + 
			len.at.res + ni.age + ni.marital.status + premium + 
			ni.gender + house.color + year + dwelling.type + coverage.type, 
			data = data.train[train,], family = binomial(link = "logit") )
			temp = step(logit, direction = "both", trace = 0, k = lambda)

		str = paste(deparse(temp$call), collapse = '')
		str = sub("~","@",str)
		str = sub(",","@",str)
		str = strsplit(str, split = "@")[[1]][2]
		str = stri_replace_all_charclass(str, "\\p{WHITE_SPACE}", "")
		str = unlist(strsplit(str, split = "[+]"))

		for(j in 1:16)
			if(names(c)[j] %in% str)
				c[j] = c[j] + 1
	}
	return(c)
}

################################################
# step 2: greedy selection of two way interactions
# factor(claim.ind), sales.channel, credit, factor(zip.code), ni.marital.status 
# n.adults, n.children, tenure, len.at.res, ni.age

V0 = c("factor(claim.ind):sales.channel", "factor(claim.ind):credit", 
		"factor(claim.ind):factor(zip.code)", 
		"factor(claim.ind):ni.marital.status",
		"sales.channel:credit", "sales.channel:factor(zip.code)", 
		"sales.channel:ni.marital.status",
		"credit:factor(zip.code)", "credit:ni.marital.status", 
		"factor(zip.code):ni.marital.status",
		"n.adults:n.children", "n.adults:tenure", 
		"n.adults:len.at.res", "n.adults:ni.age",
		"n.children:tenure","n.children:len.at.res","n.children:ni.age",
		"tenure:len.at.res", "tenure:ni.age", "len.at.res:ni.age")

M0.forluma = paste("n.adults+n.children+tenure+sales.channel+credit+factor(zip.code)", "+len.at.res+ni.age+factor(claim.ind)+ni.marital.status")

# this function returns a formula of the best model with two way interactions

InteractSelection = function(V0, M0.forluma, M0.err = 0.75, M0.auc = 0.5, 
												cv.k = 500, cv.ratio = 0.2){

	express = c("glm(cancel ~", ",data = data.train[train,], family = binomial)")

	while(TRUE){

		R = array(NA, dim = c(length(V0), 2))

		for(i in 1:length(V0)){
			n.data = nrow(data.train)

			err = numeric(cv.k)
			cstat = numeric(cv.k)

			for (j in 1:cv.k){
				train = sample(1:n.data)[1:round(n.data*cv.ratio)]
				M1 = eval(parse(paste(express[1], M0.forluma, "+", V0[i], express[2])))
				pred.p = predict(M1, data.train[-train,], type = "response")

				input = prediction(pred.p, data.train[-train,]$cancel)
				auc = performance(input, "auc")
				cstat[j] = unlist(auc@y.values) 

				pred.c = ifelse(pred.p > 1/6,1,0)
				err[j] = with(data.train[-train,], misclassification(cancel,pred.c))
			}

			M1.err = mean(err)
			M1.auc = mean(cstat)

			if(M1.err < M0.err & M1.auc > M0.auc) {
				R[i,1] = M1.err
				R[i,2] = M1.auc
			}

		}

		if( length(na.omit(R[,1])) > 0 ){
			idx = which.min(R[,1] - R[,2])
			M0.forluma = paste(M0.forluma, "+", V0[idx])
			M0.err = R[idx,1]
			M0.auc = R[idx,2]
			V0 = V0[-idx]
		}
		else break
	}
	return(M0.forluma)
}
