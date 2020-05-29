
#' dimsumms__replicate_bottleneck
#'
#' Simulate a replicate bottleneck.
#'
#' @param input_file path to input file (required)
#' @param outpath output path for plots and saved objects (required)
#' @param Bottleneck.Factor bottleneck factor (default:0.01)
#' @param Sequencing.Error sequencing error (default:0.02)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms__replicate_bottleneck <- function(
  input_file,
  outpath,
  Bottleneck.Factor=0.01,
  Sequencing.Error=0.02
  ){

  # load("/users/project/prj004631/afaure/DMS/dimsumrun_BB_TARDBP_290/BB_TARDBP_290_2017-06-13/BB_TARDBP_290_2017-06-13_variant_data_merge.RData")

  Bottlenecked.Data <- as.data.frame(fread(input_file))

  set.seed(123)
  New.Input.1 <- rpois(n = length(Bottlenecked.Data$input1_e1_s0_bNA_count),
                       lambda = (Bottlenecked.Data$input1_e1_s0_bNA_count * Bottleneck.Factor))
  New.Input.1 <- sample(x = 1:nrow(Bottlenecked.Data),
                        size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count),
                        prob = New.Input.1,
                        replace = TRUE)
  New.Input.1 <- factor(New.Input.1, levels = 1:nrow(Bottlenecked.Data))
  New.Input.1 <- as.numeric(table(New.Input.1))

  New.Input.2 <- rpois(n = length(Bottlenecked.Data$input2_e2_s0_bNA_count),
                       lambda = (Bottlenecked.Data$input2_e2_s0_bNA_count * Bottleneck.Factor))
  New.Input.2 <- sample(x = 1:nrow(Bottlenecked.Data),
                        size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count),
                        prob = New.Input.2,
                        replace = TRUE)
  New.Input.2 <- factor(New.Input.2, levels = 1:nrow(Bottlenecked.Data))
  New.Input.2 <- as.numeric(table(New.Input.2))

  New.Input.3 <- rpois(n = length(Bottlenecked.Data$input3_e3_s0_bNA_count),
                       lambda = (Bottlenecked.Data$input3_e3_s0_bNA_count * Bottleneck.Factor))
  New.Input.3 <- sample(x = 1:nrow(Bottlenecked.Data),
                        size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count),
                        prob = New.Input.3,
                        replace = TRUE)
  New.Input.3 <- factor(New.Input.3, levels = 1:nrow(Bottlenecked.Data))
  New.Input.3 <- as.numeric(table(New.Input.3))

  New.Input.4 <- rpois(n = length(Bottlenecked.Data$input4_e4_s0_bNA_count),
                       lambda = (Bottlenecked.Data$input4_e4_s0_bNA_count * Bottleneck.Factor))
  New.Input.4 <- sample(x = 1:nrow(Bottlenecked.Data),
                        size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count),
                        prob = New.Input.4,
                        replace = TRUE)
  New.Input.4 <- factor(New.Input.4, levels = 1:nrow(Bottlenecked.Data))
  New.Input.4 <- as.numeric(table(New.Input.4))

  New.Output.1 <- Bottlenecked.Data$output1A_e1_s1_b1_count
  New.Output.1[which(New.Input.1 == 0)] <- 0

  New.Output.1 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count),
                         prob = New.Output.1,
                         replace = TRUE)
  New.Output.1 <- factor(New.Output.1, levels = 1:nrow(Bottlenecked.Data))
  New.Output.1 <- as.numeric(table(New.Output.1))

  New.Output.2 <- Bottlenecked.Data$output2A_e2_s1_b1_count
  New.Output.2[which(New.Input.2 == 0)] <- 0

  New.Output.2 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count),
                         prob = New.Output.2,
                         replace = TRUE)
  New.Output.2 <- factor(New.Output.2, levels = 1:nrow(Bottlenecked.Data))
  New.Output.2 <- as.numeric(table(New.Output.2))

  New.Output.3 <- Bottlenecked.Data$output3A_e3_s1_b1_count
  New.Output.3[which(New.Input.3 == 0)] <- 0

  New.Output.3 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count),
                         prob = New.Output.3,
                         replace = TRUE)
  New.Output.3 <- factor(New.Output.3, levels = 1:nrow(Bottlenecked.Data))
  New.Output.3 <- as.numeric(table(New.Output.3))


  New.Output.4 <- Bottlenecked.Data$output4A_e4_s1_b1_count
  New.Output.4[which(New.Input.4 == 0)] <- 0

  New.Output.4 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count),
                         prob = New.Output.4,
                         replace = TRUE)
  New.Output.4 <- factor(New.Output.4, levels = 1:nrow(Bottlenecked.Data))
  New.Output.4 <- as.numeric(table(New.Output.4))

  Bottlenecked.Data$input1_e1_s0_bNA_count <- New.Input.1
  Bottlenecked.Data$input2_e2_s0_bNA_count <- New.Input.2
  Bottlenecked.Data$input3_e3_s0_bNA_count <- New.Input.3
  Bottlenecked.Data$input4_e4_s0_bNA_count <- New.Input.4

  Bottlenecked.Data$output1A_e1_s1_b1_count <- New.Output.1
  Bottlenecked.Data$output2A_e2_s1_b1_count <- New.Output.2
  Bottlenecked.Data$output3A_e3_s1_b1_count <- New.Output.3
  Bottlenecked.Data$output4A_e4_s1_b1_count <- New.Output.4

  New.Output.1 <- Bottlenecked.Data$output1B_e1_s1_b2_count
  New.Output.1[which(New.Input.1 == 0)] <- 0

  New.Output.1 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count),
                         prob = New.Output.1,
                         replace = TRUE)
  New.Output.1 <- factor(New.Output.1, levels = 1:nrow(Bottlenecked.Data))
  New.Output.1 <- as.numeric(table(New.Output.1))

  New.Output.2 <- Bottlenecked.Data$output2B_e2_s1_b2_count
  New.Output.2[which(New.Input.2 == 0)] <- 0

  New.Output.2 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count),
                         prob = New.Output.2,
                         replace = TRUE)
  New.Output.2 <- factor(New.Output.2, levels = 1:nrow(Bottlenecked.Data))
  New.Output.2 <- as.numeric(table(New.Output.2))

  New.Output.3 <- Bottlenecked.Data$output3B_e3_s1_b2_count
  New.Output.3[which(New.Input.3 == 0)] <- 0

  New.Output.3 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count),
                         prob = New.Output.3,
                         replace = TRUE)
  New.Output.3 <- factor(New.Output.3, levels = 1:nrow(Bottlenecked.Data))
  New.Output.3 <- as.numeric(table(New.Output.3))


  New.Output.4 <- Bottlenecked.Data$output4B_e4_s1_b2_count
  New.Output.4[which(New.Input.4 == 0)] <- 0

  New.Output.4 <- sample(x = 1:nrow(Bottlenecked.Data),
                         size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count),
                         prob = New.Output.4,
                         replace = TRUE)
  New.Output.4 <- factor(New.Output.4, levels = 1:nrow(Bottlenecked.Data))
  New.Output.4 <- as.numeric(table(New.Output.4))

  Bottlenecked.Data$output1B_e1_s1_b2_count <- New.Output.1
  Bottlenecked.Data$output2B_e2_s1_b2_count <- New.Output.2
  Bottlenecked.Data$output3B_e3_s1_b2_count <- New.Output.3
  Bottlenecked.Data$output4B_e4_s1_b2_count <- New.Output.4

  Idx.WT <- which(Bottlenecked.Data$Nham_nt == 0)
  Idx.Singles <- which(Bottlenecked.Data$Nham_nt == 1)
  Idx.Doubles <- which(Bottlenecked.Data$Nham_nt == 2)
  Idx.Triples <- which(Bottlenecked.Data$Nham_nt == 3)
  Idx.Quadruples <- which(Bottlenecked.Data$Nham_nt == 4)

  # sequencing errors input 1

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Quadruples] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.WT] <- Bottlenecked.Data$input1_e1_s0_bNA_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors input 2

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Quadruples] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.WT] <- Bottlenecked.Data$input2_e2_s0_bNA_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors input 3

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Quadruples] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.WT] <- Bottlenecked.Data$input3_e3_s0_bNA_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors input 4

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Quadruples] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.WT] <- Bottlenecked.Data$input4_e4_s0_bNA_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 1

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Quadruples] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))))
  Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.WT] <- Bottlenecked.Data$output1A_e1_s1_b1_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 2

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Quadruples] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.WT] <- Bottlenecked.Data$output2A_e2_s1_b1_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 3

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Quadruples] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.WT] <- Bottlenecked.Data$output3A_e3_s1_b1_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 4

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Quadruples] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples] - Counts.To.Subtract.From.Triples

  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles

  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.WT] <- Bottlenecked.Data$output4A_e4_s1_b1_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 1

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Quadruples] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples] - Counts.To.Subtract.From.Triples


  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles


  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))))
  Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.WT] <- Bottlenecked.Data$output1B_e1_s1_b2_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 2

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Quadruples] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples] - Counts.To.Subtract.From.Triples


  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles


  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.WT] <- Bottlenecked.Data$output2B_e2_s1_b2_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 3

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Quadruples] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples] - Counts.To.Subtract.From.Triples


  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles


  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.WT] <- Bottlenecked.Data$output3B_e3_s1_b2_count[Idx.WT] - Counts.To.Subtract.From.WT

  # sequencing errors output 4

  Counts.To.Add.To.Quadruples <- as.vector(rmultinom(n = 1,
                                                     size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                     prob = rep(1, times = length(Idx.Quadruples))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Quadruples] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Quadruples] + Counts.To.Add.To.Quadruples


  Counts.To.Subtract.From.Triples <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples] - Counts.To.Subtract.From.Triples


  Counts.To.Add.To.Triples <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Triples))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Triples] + Counts.To.Add.To.Triples


  Counts.To.Subtract.From.Doubles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles] - Counts.To.Subtract.From.Doubles


  Counts.To.Add.To.Doubles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Doubles))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Doubles] + Counts.To.Add.To.Doubles


  Counts.To.Subtract.From.Singles <- as.vector(rmultinom(n = 1,
                                                         size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles]) * Sequencing.Error,
                                                         prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles] - Counts.To.Subtract.From.Singles


  Counts.To.Add.To.Singles <- as.vector(rmultinom(n = 1,
                                                  size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                  prob = rep(1, times = length(Idx.Singles))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.Singles] + Counts.To.Add.To.Singles


  Counts.To.Subtract.From.WT <- as.vector(rmultinom(n = 1,
                                                    size = sum(Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.WT]) * Sequencing.Error,
                                                    prob = rep(1, times = length(Idx.WT))
  ))
  Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.WT] <- Bottlenecked.Data$output4B_e4_s1_b2_count[Idx.WT] - Counts.To.Subtract.From.WT

  #Finalise
  variant_data_merge <- as.data.table(Bottlenecked.Data)
  #Replace negative counts with 0
  for(j in grep("_count$", names(variant_data_merge))){
    set(variant_data_merge, i=which(variant_data_merge[[j]]<0), j=j, value=0)
  }
  #Remove variants with zero counts in all samples
  variant_data_merge <- variant_data_merge[apply(variant_data_merge[,.SD,.SDcols = names(variant_data_merge)[grep("_count$", names(variant_data_merge))]], 1, sum)!=0,]
  #Reformat
  variant_data_merge <- variant_data_merge[,.SD,,.SDcols = grep("nt_seq|_count$", names(variant_data_merge))]
  names(variant_data_merge) <- sapply(strsplit(names(variant_data_merge), "_e"), '[', 1)

  #Save
  OutFileName <- file.path(outpath, paste0("BB_TARDBP_290_2017-06-13_ReplicateBottleneck_", substr(as.character(Bottleneck.Factor), 3, 4), "_t0_variant_data_merge.tsv"))
  write.table(variant_data_merge, OutFileName, sep = "\t", quote = F, row.names = F)

}

