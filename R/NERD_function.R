#' The Nosocomial Evolution of Resistance Detector (NERD)
#'
#' The NERD package analyses temporal trends in antibiotic resistance data.
#' @param
#' times The date information. Should be stored in non-negative integers (e.g. days).
#' @param
#' res_info The resistance information. Should be stored on a discrete scale from 1 to some value m. Values below 0 are discarded. For us: load "CData" which stores evolution information and select appropriate column to extract antibiotic of interest.
#' @param
#' param The parameter for the Poisson part of the test. Defaults to 365 days = 1 year.
#' @keywords Antibiotic Resistance
#' @export
#' @examples
#' NERD_function()

NERD_function <- function(times,res_info,param=365) { 
	
	N <- length(times);
	NN <- length(res_info);
	
	if (N != NN) {
		stop("Time information and resistance information have the same length")
	}
	
	m <- max(res_info);
	m_counts <- rep(0,m);
	for (i in 1:m) {
		m_counts[i] <- sum(res_info == i);
	}
	min_m <- min(m_counts);
	if (min_m < 1) {
		stop("Must be a postive number of counts in each resistance category")
	}
		
	
	

	
	times <- times[res_info > 0];
	res_info <- res_info[res_info > 0];
	ptimes1 <- array( rep(times, N), c(N,N));
	ptimes2 <- t(ptimes1);
	ptimes <- ptimes1 - ptimes2;
	ptimes[ptimes < 0] <- 0;
	#m = max(res_info);
	
	Markov = array( 0,c(m,m))
	for (i in 1:m) {
		for (j in 1:m) {
			cols <- (res_info == j);
			rows <- (res_info == i);
			pt <- ptimes[rows,cols]; # raw pairwise time differences
			ptw <- dpois(pt, param, log = FALSE); # weighted by Poisson
			Markov[j,i] <- sum(ptw);
		}
	}
	data_sums = rep(0,m);
	nums = rep(0,m);
	for (i in 1:m) {
		data_sums[i] = sum(Markov[i,]);
		nums[i] = sum(res_info == i);
	}
	Counts = nums/sum(res_info >0); # gets proportion of each status
	markov_data <- array(0,c(m,m));
	for (i in 1:m) {
		markov_data[i,] <- Markov[i,]/data_sums[i];
	}
	Counts1 <- t(array( rep(Counts), c(m,m)));
	chi <- (markov_data - Counts1)^2;
	chi <- chi/Counts1;
	chi_data <- sum(chi)







	# then we run simulations to compare

	M <- 1000; # the number of randomized strings we compare to
	chis <- rep(0, M);
	for (jj in 1:M) {
		
		res_info <- rep(1,nums[1]);
		for (i in 2:m) {
			res_info <- c(res_info, rep(i,nums[i]));
		}
		res_info <- sample(res_info)
		ptimes1 <- array( rep(times, N), c(N,N));
		ptimes2 <- t(ptimes1);
		ptimes <- ptimes1 - ptimes2;
		ptimes[ptimes < 0] <- 0
		Markov = array( 0,c(m,m))
		for (i in 1:m) {
			for (j in 1:m) {
				cols <- (res_info == j);
				rows <- (res_info == i);
				pt <- ptimes[rows,cols];
				ptw <- dpois(pt, param, log = FALSE)
				Markov[j,i] <- sum(ptw);
			}
		}

		sums = rep(0,m);
		counts = rep(0,m);
		for (i in 1:m) {
			sums[i] = sum(Markov[i,]);
			counts[i] = sum(res_info == i);
		}
		counts = counts/sum(res_info > 0);
		Markov1 = array(0,c(m,m))
		for (i in 1:m) {
			Markov1[i,] <- Markov[i,]/sums[i];
		}
		Counts1 <- t(array( rep(Counts), c(3,3)));
		chi <- (Markov1 - Counts1)^2;
		chi <- chi/Counts1;
		chi <- sum(chi)
		chis[jj] <- chi;
	}

	chis1 <- sort(chis);
	pvalue <- sum(chis > chi_data)/M;

	markov_data;
	Counts;
	pvalue;
}



