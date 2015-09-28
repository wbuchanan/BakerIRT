package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * A program to estimate an examinee's ability under
 * the one-, two-, or three-parameter logistic item
 * characteristic curve model.
 */
public class MLEPersonParameter {



		static final int NITEM = 10;
		static final double BIGT = 0.5;
		public MLEPersonParameter(String[] args) {
			double[] cpt = new double[NITEM];
			cpt[0] =  0.489; // read canned data
			cpt[1] = -0.869;
			cpt[2] = -1.247;
			cpt[3] =  0.595;
			cpt[4] =  0.126;
			cpt[5] =  0.469;
			cpt[6] =  0.058;
			cpt[7] =  0.665;
			cpt[8] = -1.292;
			cpt[9] = -0.490;
			double[] a = { 0.997, 0.874, 1.110, 1.435, 1.351,
					1.501, 1.955, 0.906, 1.897, 1.146 };
			double[] c = { .15, .21, .23, .01, .21,
					.02, .27, .05, .31, .27 };
			double[] uij = { 1.0, 1.0, 1.0, 0.0, 0.0,
					0.0, 1.0, 0.0, 1.0, 1.0 };
			int trall = -1; // establish traces
			int trmle = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trmle = 0;
			}
			else{
				String yorn = JOptionPane.showInputDialog
						("Trace MLE? Y/N");
				if (yorn.equalsIgnoreCase("y")) trmle = 0;
			} // end if (yn ... else
			String smaxit = JOptionPane.showInputDialog
					("Enter maximum number of iterations to do. ");
			int maxit = Integer.parseInt(smaxit);
			String snpara = JOptionPane.showInputDialog
					("Enter number of parameters in ICC model. ");
			int npara = Integer.parseInt(snpara);
			System.out.println("All= " + trall + "  MLE= " + trmle +
					"  Max Iterations= " + maxit + "  ICC Model= " + npara);
			if (npara == 1) {
				for (int i = 0; i <= NITEM - 1; i++) {
					a[i] = 1.0;
				} // end for
			}
			else {
				// normal processing
			} // end if (npara ... else
			double theta = 0.0;
			// estimate ability
			// ability estimation subroutine
			int lastnit = 1;
			for (int nit = 1; nit <= maxit; nit++) {
				lastnit = nit;
				double sumnum = 0.0;
				double sumdem = 0.0; // clear sums
				System.out.println("");
				System.out.println("Iteration= " + nit);
				for (int i = 0; i <= NITEM - 1; i++) { // item loop
					// calculate phat
					double dev = cpt[i] + a[i] * theta;
					double phat = 1.0 / (1.0 + Math.exp(-dev));
					if (npara != 3) { // one- and two-parameter ICC models
						double wij = phat * (1.0 - phat);
						double vij = uij[i] - phat;
						sumnum += a[i] * vij;
						sumdem += a[i] * a[i] * wij;
						if (trall == 0) {
							System.out.println("");
							System.out.println("Item= " + (i + 1));
							System.out.println("Intercept= " + cpt[i] +
									"  Slope= " + a[i] + "  Theta= " + theta);
							System.out.println("Uij= " + uij[i] +
									"  Dev= " + dev + "  Phat= " + phat);
							System.out.println("Wij= " + wij + "  Vij= " + vij);
							System.out.println("Numerator Sum= " + sumnum +
									"  Denominator Sum= " + sumdem);
						}
						else {
							// normal processing
						} // end if (trall ... else
					}
					else { // three-parameter ICC model
						double pt = c[i] + (1.0 - c[i]) * phat;
						// protect against divide by zero
						if (pt < 0.00001) pt = 0.00001;
						if (pt > 0.99999) pt = 0.99999;
						double wij = pt * (1.0 - pt);
						double vij = uij[i] - pt;
						double psp = phat / pt;
						sumnum += a[i] * vij * psp;
						sumdem += a[i] * a[i] * wij * psp * psp;
						if (trall == 0) {
							System.out.println("");
							System.out.println("Item= " + (i + 1) +
									"  C= " + c[i]);
							System.out.println("Intercept= " + cpt[i] +
									"  Slope= " + a[i] + "  Theta= " + theta);
							System.out.println("Uij= " + uij[i] +
									"  Dev= " + dev + "  P*= " + phat);
							System.out.println("P= " + pt + "  P*/P= " + psp);
							System.out.println("Wij= " + wij + "  Vij= " + vij);
							System.out.println("Numerator Sum= " + sumnum +
									"  Denominator Sum= " + sumdem);
						}
						else {
							// normal processing
						} // end if (trall ... else
					} // end if (npara ... else
				} // end for (int i ... (BASIC 2370)
				double delta = sumnum / sumdem;
				if (trall == 0 || trmle == 0) {
					System.out.println("");
					System.out.println("Change in Theta= " + delta);
				}
				else {
					// normal processing
				} // end if (trall ... else
				if (Math.abs(delta) > BIGT) {
					if (delta > 0.0) {
						delta = BIGT;
					}
					else {
						delta = - BIGT;
					} // end if (delta ... else
				}
				else {
					// normal processing
				} // end if (Math ... else
				theta += delta;
				if (trall == 0 || trmle == 0) {
					System.out.println("");
					System.out.println("Theta= " + theta +
							"  Change= " + delta);
				}
				else {
					// normal processing
				} // end if
				if (Math.abs(delta) < 0.05) {
					break; // exit the for statement
				}
				else {
					// normal processing
				} // end if (Math ... else
			} // end for (int nit ...
			if (lastnit == maxit) System.out.println
					("Reached max iteration. ");
			System.out.println("");
			System.out.println("Estimated Ability= " + theta);
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
