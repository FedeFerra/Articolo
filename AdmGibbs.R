Adm.Gibbs <- function(dtm, K, iteration, alpha = NULL, beta = NULL, verbose = TRUE, post = TRUE){

  N_d_w <- as.matrix(dtm)
  V <- ncol(dtm)
  D <- nrow(dtm)

  alpha <- 1/K 
  beta <- 1/V
  z <- rep(0,D)
  #begin
  m_z <- rep(0,K)
  n_z <- rep(0,K)
  n_w_z <- matrix(0,nrow = V,ncol=K)
  N_d <- as.vector(rowSums(N_d_w))

  for(d in 1:D){
    z[d] <- which(rmultinom(1, 1, rep(1/K, K)) == 1)
    m_z[z[d]] <- m_z[z[d]] + 1
    n_z[z[d]] <- n_z[z[d]] + N_d[d]
    n_w_z[, z[d]] <- n_w_z[, z[d]] + N_d_w[d, ]
  }

  
  occurrences <- list()
  doc.topic <- list()
  posterior <- list()
  kernel.matrix <- matrix(NA, D, K)

  for(i in 1:iteration){
    #Andranno salvate gli output per ogni generazione
    for(d in 1:D){
      m_z[z[d]] <- m_z[z[d]] - 1
      n_z[z[d]] <- n_z[z[d]] - N_d[d]
      n_w_z[, z[d]] <- n_w_z[, z[d]] - N_d_w[d, ]
      #Estrarre dalla full conditional
      index <- which(N_d_w[d, ] != 0)
      vector.constant <- (m_z + rep(alpha, length(m_z)))/rep(D - 1 + K*alpha, length(m_z))
    
      conteggi <- N_d_w[d, index]
      J <- sequence(conteggi)
    
      matrix.J <- matrix(J, nrow = length(J), ncol = K)
      matrix.Beta <- matrix(beta, nrow = length(J), ncol = K)
      matrix.One <- matrix(1, nrow = length(J), ncol = K)
    
      #Ho aggiunto un as.matrix qui
      matrix.count <- apply(as.matrix(n_w_z[index, ]), 2, function(x) rep(x, as.vector(conteggi)))
      #La riga qui sotto serve per risolvere il problema nel caso in cui ci sia n_w_z 1xK
      matrix.count <- matrix(as.vector(matrix.count), nrow = sum(conteggi), ncol = K)

      matrix.sum <- matrix.J + matrix.Beta + matrix.count - matrix.One
      vector.numerator <- apply(matrix.sum, 2, function(x) prod(x))
    
      I <- sequence(N_d[d])
    
      matrix.I <- matrix(I, nrow = length(I), ncol = K)
      matrix.Vbeta <- matrix(V*beta, nrow = length(I), ncol = K)
      matrix.One <- matrix(1, nrow = length(I), ncol = K)
      matrix.count <- matrix(n_z, nrow = length(I), ncol = K, byrow = TRUE)
    
      matrix.sum <- matrix.count + matrix.Vbeta + matrix.I - matrix.One
      vector.denominator <- apply(matrix.sum, 2, function(x) prod(x))
    
      kernel <- vector.constant*(vector.numerator/vector.denominator)
      if(post){
      kernel.matrix[d, ] <- kernel
      }
      
      z[d] <- sample(1:K, 1, prob = kernel)
      #Aumento i conteggi per topic e e parole totali
      m_z[z[d]] <- m_z[z[d]]+1
      n_z[z[d]] <- n_z[z[d]]+N_d[d]
      }
      #Aumento i conteggi per parole nel documento corrente
      occurrences[[i]] <- n_w_z
      doc.topic[[i]] <- z
      if(post){
        posterior[[i]] <- kernel.matrix
      }
      if(verbose & i%%500 == 0){
        cat(i)
      } 
    }
list(occurrences.matrix = occurrences, topic.assignment = doc.topic, posterior.estimates = posterior)
}

