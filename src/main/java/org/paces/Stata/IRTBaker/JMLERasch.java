package org.paces.Stata.IRTBaker;

import javax.swing.*;


/**
 * A program to implement the Birnbaum joint maximum likelihood
 * estimation paradigm for the Rasch model.
 */
public class JMLERasch {



		static final int NITEM = 10;
		static final int MAXSCORE = 9;
		static final int MAXCYCLE = 10;
		static final int MAXIT = 5;

		public JMLERasch(String[] args) {

			// Raw score frequencies
			double[] fdg = new double[MAXSCORE]; // read canned data, ...
			fdg[0] =  53.0; // frequency of raw scores for 1 to MAXSCORE
			fdg[1] =  85.0;
			fdg[2] = 127.0;
			fdg[3] = 148.0;
			fdg[4] = 163.0;
			fdg[5] = 120.0;
			fdg[6] = 135.0;
			fdg[7] =  95.0;
			fdg[8] =  50.0;


			// Raw Scores
			double[] s = { 497.0, 622.0, 366.0, 522.0, 508.0,
					381.0, 436.0, 526.0, 400.0, 628.0 };

			// Number of items in the model
			double[] b = new double[NITEM];

			// Unique total scores
			double[] theta = new double[MAXSCORE];



			int trall = -1;
			int tritm = -1;
			int trabl = -1;

			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				tritm = 0;
				trabl = 0;
			}
			else{
				String ynitem = JOptionPane.showInputDialog
						("Trace item MLE? Y/N");
				if (ynitem.equalsIgnoreCase("y")) {
					tritm = 0;
				}
				String ynability = JOptionPane.showInputDialog
						("Trace ability MLE? Y/N");
				if (ynability.equalsIgnoreCase("y")) {
					trabl = 0;
				}
			} // end if (yn ... else


			// Iterate over the number of unique total score values and sum
			// the cell frequencies.  Return the scalar containing the number
			// of observations.
			double n = 0.0;
			for (int j = 0; j <= MAXSCORE - 1; j++) {
				n += fdg[j];
			}

			if (trall == 0) {
				System.out.println("N= " + n);
			}

			// Iterate over items and
			double sumi = 0.0;
			for (int i = 0; i <= NITEM - 1; i++) {

				b[i] = Math.log((n - s[i]) / s[i]);
				sumi += b[i];
				if (trall == 0) {
					System.out.println("SumI= " + sumi + "  B(" + (i + 1) +
							")= " + b[i]);
				}
			} // end for (int i ...

			double meanl = sumi / NITEM; // casting
			if (trall == 0) {
				System.out.println("Mean of log terms= " + meanl);
			}


			for (int i = 0; i <= NITEM - 1; i++) {
				b[i] -= meanl;
				if (trall == 0) {
					System.out.println("Initial B(" + (i + 1) + ")= " + b[i]);
				}
			}

			/* Preliminary estimate of Theta is
			 */
			for (int j = 0; j <= MAXSCORE - 1; j++) {
				theta[j] = Math.log((j + 1.0) / (NITEM - (j + 1.0)));
				// casting
				if (trall == 0) {
					System.out.println("Initial theta(" + (j + 1) + ")= " +
							theta[j]);
				}
			} // BASIC 3260


			double bbold = 9999.99; // initialize bbold

			// Loop over iteration cycles
			for (int k = 1; k <= MAXCYCLE; k++){
				System.out.println("JMLE cycle= " + k);

				// Loop over items
				double sumb = 0.0;
				for (int i = 0; i <= NITEM - 1; i++) { // BASIC 1000

					for (int kk = 1; kk <= MAXIT; kk++) {
						double sumfdgp = 0.0;
						double sumfdgpq = 0.0;

						for (int j = 0; j <= MAXSCORE - 1; j++) {
							if (trall == 0) {
								System.out.println("Item= " + i +
										"  NR Iteration= " + kk +
										"  Raw Score= " + (j + 1));
							}

							double dev = theta[j] - b[i];
							double pij = 1.0 / (1 + Math.exp(-dev));
							sumfdgp += fdg[j] * pij;
							sumfdgpq += fdg[j] * pij * (1.0 - pij);
							if (trall == 0) {
								System.out.println("Dev= " + dev +
										"  Pij= " + pij +
										"  SumFDGP= " + sumfdgp +
										"  SumFDGPQ= " + sumfdgpq);
							}
							else {
								// normal processing
							} // end if (trall ... else
						} // end for (int j ... (BASIC 1120)
						double deltab = (s[i] - sumfdgp) / sumfdgpq;
						b[i] = b[i] - deltab;
						if (trall == 0 || tritm == 0) {
							System.out.println("DeltaB= " + deltab +
									"  B(" + (i + 1) + ")= " + b[i]);
						}
						else {
							// normal processing
						} // end if (trall ... else
						if (Math.abs(deltab) < 0.05) {
							break; // exit the for statement
						}
						else {
							// normal processing
						} // end if (Math ... else
					} // end for (int kk ... (BASIC 1170)
					sumb += b[i];
				} // end for (int i ... (BASIC 1190)
				double bbar = sumb / NITEM; // casting
				for (int i = 0; i <= NITEM - 1; i++) { // correction
					b[i] -= bbar;
				}
				if (trall == 0 || tritm == 0) {
					System.out.println("SumB= " + sumb +
							"  BBar= " + bbar +
							"  Nitem= " + NITEM);
				}
				else {
					// normal processing
				} // end if (trall ... else
				for (int i = 0; i <= NITEM - 1; i++) {
					System.out.println("B(" + (i + 1) + ")= " + b[i]);
				}
				if (Math.abs(bbar - bbold) > 0.05) {
					// Ability Estimation by Raw Score
					for (int j = 0; j <= MAXSCORE - 1; j++) { // BASIC 2000
						for (int kk = 1; kk <= MAXIT; kk++) {
							double sump = 0.0;
							double sumpq = 0.0;
							for (int i = 0; i <= NITEM - 1; i++) {
								double tdev = theta[j] - b[i];
								double tpij = 1.0 / (1.0 + Math.exp(-tdev));
								sump += tpij;
								sumpq += tpij * (1.0 - tpij);
								if (trall == 0) {
									System.out.println("Raw Score= " + (j + 1) +
											"  NR Iteration= " + kk +
											"  Item= " + (i + 1));
									System.out.println("Deve= " + tdev +
											"  Pij= " + tpij +
											"  SumP= " + sump +
											"  SumPQ= " + sumpq);
								}
							} // end for (int i ... (BASIC 2110)
							double delta = (j + 1.0 - sump) / sumpq; // casting
							theta[j] += delta;
							if (trall == 0 || trabl == 0) {
								System.out.println("Delta= " + delta +
										"  Theta(" + (j + 1) + ")= " + theta[j]);
							} // end if (trall ...
							if (Math.abs(delta) < 0.05) {
								break; // exit the for (int kk ...
							}
						} // end for (int kk ... (BASIC 2160)
						if (trall == 0 || trabl == 0) {
							System.out.println("Theta(" + (j + 1) + ")= " + theta[j]);
						} // end if (trall ...
					} // end for (int j ... (BASIC 2180)
				}
				else {
					System.out.println("Reached convergence. ");
					break; // exit the for statement
				} // end if (Math.abs( ... else
				if (k == MAXCYCLE) {
					System.out.println("Max cycles reached. ");
				}
			} //end for (int k ... (BASIC 390)
			// correct item difficulty for bias
			System.out.println("Item Difficulty Corrected for Bias");
			for (int i = 0; i <= NITEM -1; i++) {
				b[i] = b[i] * ((NITEM - 1.0) / NITEM); // casting
				System.out.println("B(" + (i + 1) + ")= " + b[i]);
			}
			// Ability Estimation by Raw Score
			for (int j = 0; j <= MAXSCORE - 1; j++) { // BASIC 2000
				for (int kk = 1; kk <= MAXIT; kk++) {
					double sump = 0.0;
					double sumpq = 0.0;
					for (int i = 0; i <= NITEM - 1; i++) {
						double tdev = theta[j] - b[i];
						double tpij = 1.0 / (1.0 + Math.exp(-tdev));
						sump += tpij;
						sumpq += tpij * (1.0 - tpij);
						if (trall == 0) {
							System.out.println("Raw Score= " + (j + 1) +
									"  NR Iteration= " + kk +
									"  Item= " + (i + 1));
							System.out.println("Deve= " + tdev +
									"  Pij= " + tpij +
									"  SumP= " + sump +
									"  SumPQ= " + sumpq);
						}
					} // end for (int i ... (BASIC 2110)
					double delta = (j + 1.0 - sump) / sumpq; // casting
					theta[j] += delta;
					if (trall == 0 || trabl == 0) {
						System.out.println("Delta= " + delta +
								"  Theta(" + (j + 1) + ")= " + theta[j]);
					} // end if (trall ...
					if (Math.abs(delta) < 0.05) {
						break; // exit the for (int kk ...
					}
				} // end for (int kk ... (BASIC 2160)
				if (trall == 0 || trabl == 0) {
					System.out.println("Theta(" + (j + 1) + ")= " + theta[j]);
				} // end if (trall ...
			} // end for (int j ... (BASIC 2180)
			System.out.println("Ability Corrected for Bias");
			for (int j = 0; j <= MAXSCORE - 1; j++) {
				theta[j] = theta[j] * (NITEM - 2.0) / NITEM; // casting
				System.out.println("Theta(" + (j + 1) + ")= " + theta[j]);
			} // end for (int ... (BASIC 370)
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
