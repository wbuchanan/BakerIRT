package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * A program to get maximum likelihood estimates of the
 * slope and intercept parameters for a single item under
 * the two-parameter logistic item characteristic curve
 * model.
 */
public class MLEItemParameters {

		public MLEItemParameters(String[] args) {
			double[] x = {-4.0000, -3.1111, -2.2222, -1.3333,
					-0.4444, 0.4444, 1.3333, 2.2222, 3.1111, 4.0000};
			double[] r = {6.0, 17.0, 20.0, 34.0, 51.0, 68.0,
					81.0, 90.0, 95.0, 97.0};
			int nxl = 9;
			double[] f = new double[10];
			for (int k = 0; k <= nxl; k++) f[k] = 100.0;
			double cpt = 0.0;
			double a = 1.0; //initial values
			int trall = -1;
			int trmle = -1; // set up trace of computations
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
			} // end if else
			String smaxit = JOptionPane.showInputDialog
					("Enter number of iterations to do. ");
			int maxit = Integer.parseInt(smaxit);
			System.out.println("All= " + trall + "  MLE= " + trmle +
					"  Max Iterations= " + maxit);
			int lastnit = 1;
			double lastsfwv2 = 0.0;
			for (int nit = 1; nit <= maxit; nit++) {
				lastnit = nit;
				System.out.println(" ");
				System.out.println("Iteration= " + nit);
				double sfw = 0.0;
				double sfwx = 0.0;
				double sfwx2 = 0.0;
				double sfwv = 0.0;
				double sfwvx = 0.0;
				double sfwv2 = 0.0;
				for (int k = 0; k <= nxl; k++) {
					double p = r[k] / f[k]; // use p for pi
					double dev = cpt + a * x[k];
					double ph = 1 / (1 + Math.exp(-dev));
					double w = ph * (1.0 - ph);
					if (w > 0.0000009) {
						double v = (p - ph) / w;
						if (trall == 0) {
							System.out.println("");
							System.out.println("x(" + (k+1) + ")= " + x[k] +
									"  p= " + p + "  p hat= " + ph);
							System.out.println("w= " + w + "  v= " + v);
						}
						else {
							// normal process
						} // end if else
						double fw = f[k] * w;
						double fwx = fw * x[k];
						double fwx2 = fwx * x[k];
						double fwv = fw * v;
						double fwvx = fwv * x[k];
						double fwv2 = fwv * v;
						sfw = sfw + fw;
						sfwx = sfwx + fwx;
						sfwx2 = sfwx2 + fwx2;
						sfwv = sfwv + fwv;
						sfwvx = sfwvx + fwvx;
						sfwv2 = sfwv2 + fwv2;
						if (trall == 0) {
							System.out.println("fw= " + fw + "  fwx= " +
									fwx);
							System.out.println("fwx2= " + fwx2 +
									"  fwv= " + fwv);
							System.out.println("fwvx= " + fwvx +
									"  fwv2= " + fwv2);
							System.out.println("sfw= " + sfw +
									"  sfwx= " + sfwx);
							System.out.println("sfwx2= " + sfwx2 +
									"  sfwv= " + sfwv);
							System.out.println("sfwvx= " + sfwvx +
									"  sfwv2= " + sfwv2);
						}
						else {
							// normal processing
						} // end if else
					}
					else {
						// normal processing
					} // end if else
				} // end for (i.e., BASIC line 1180)
				lastsfwv2 = sfwv2;
				if (sfw <= 0.0){
					System.out.println("");
					System.out.println
							("Out of bound error in iteration " + nit);
					break; // exit the for statement
				}
				else {
					// normal processing
				} // end if else
				double dm = sfw * sfwx2 - sfwx * sfwx;
				if (trall == 0) {
					System.out.println("");
					System.out.println ("Denominator= " + dm);
				}
				else {
					// normal processing
				} // end if else
				if (dm <= 0.000099) {
					System.out.println("");
					System.out.println("Small denominator in iteration "
							+ nit);
					break; // exit the for statement
				}
				else {
					// normal processing
				} // end if else
				double dcpt = (sfwv * sfwx2 - sfwvx * sfwx) / dm;
				double da = (sfw * sfwvx - sfwx * sfwv) /dm;
				cpt = cpt + dcpt;
				a = a + da;
				if (trmle == 0) {
					System.out.println("");
					System.out.println("nit= " + nit + "  Intercept= " + cpt +
							"  Change= " + dcpt);
					System.out.println("Slope= " + a + "  Change= " + da);
				}
				else {
					// normal processing
				} // end if else
				if (Math.abs(cpt) > 30.0 && Math.abs(a) > 20.0) {
					System.out.println("");
					System.out.println
							("Out of bound error in iteration " + nit);
					break; // exit the for statement
				}
				else {
					// normal processing
				} // end if else
				if (Math.abs(dcpt) <= 0.05 || Math.abs(da) <= 0.05) {
					break; // exit the for statement
				}
				else {
					// normal processing
				} // end if else
			} // end for (i.e., BASIC line 1330)
			if (lastnit == maxit) System.out.println
					("Reached maximum number of iterations");
			double diff = - cpt / a;
			int df = nxl - 2;
			System.out.println("");
			System.out.println("Intercept= " + cpt +
					"  Slope= " + a);
			System.out.println("Difficulty= " + diff +
					"  Discrimination= " + a);
			System.out.println("Chi-square= " + lastsfwv2 +
					"  df= " + df);
			System.exit(0); // last line
		} // end public static
	} // end public class
