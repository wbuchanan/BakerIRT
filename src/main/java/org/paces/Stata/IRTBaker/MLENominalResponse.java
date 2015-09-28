package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class MLENominalResponse {


/**
 * A program to estimate item parameters under nominal response
 * model for a three-category item with seven ability groups.
 */

		static final int MRC = 3; // nuMber of Response Categories
		static final int NTP = 2; // Number of Transformed Parameters
		static final int MAXIT = 5; // MAXimum number of ITerations
		static final int MAXGPS = 7; // number of theta Points on Scale
		public MLENominalResponse(String[] args) {
			this.nominalResponseItemParameters();
			this.nominalResponsePersonParameter();
		}

		public void nominalResponseItemParameters() {
			double[] p = new double[MRC];
			double[] cp = new double[NTP];
			double[] ap = new double[NTP];
			double[] z = new double[MRC];
			double[] ez = new double[MRC];
			double[] t = new double[NTP*2];
			double[] term = new double[NTP*2];
			double[] fds = new double[NTP*2];
			double[] pdelta = new double[NTP*2];
			double[][] tm = new double[NTP][NTP];
			double[][] w = new double[NTP][NTP];
			double[][] mtrx = new double[NTP*2][NTP*2];
			// read canned data (BASIC 500)
			double[] cpt = new double[MRC];
			cpt[0] = -0.80;
			cpt[1] = 0.35;
			cpt[2] = 0.45;
			double[] a = new double[MRC];
			a[0] = -0.30;
			a[1] = -0.60;
			a[2] = 0.90;
			double[] theta = new double[MAXGPS];
			theta[0] = -3.0;
			theta[1] = -2.0;
			theta[2] = -1.0;
			theta[3] = 0.0;
			theta[4] = 1.0;
			theta[5] = 2.0;
			theta[6] = 3.0;
			double[] f = { 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
					1000.0, 1000.0 };
			double[][] r = {
					{ 57.0, 938.0, 5.0 },
					{ 97.0, 873.0, 30.0 },
					{ 143.0, 707.0, 150.0 },
					{ 139.0, 377.0, 484.0 },
					{ 71.0, 106.0, 823.0 },
					{ 25.0, 20.0, 955.0 },
					{ 7.0, 4.0, 989.0 },
			};
			int trall = -1; // set traces
			int trmle = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trmle = 0;
			}
			else{
				String ynmle = JOptionPane.showInputDialog
						("Trace MLE terms? Y/N");
				if (ynmle.equalsIgnoreCase("y")) {
					trmle = 0;
				}
				else {
					// normal processing
				} // end if (ynmle ...
			} // end if (yn ... else
			System.out.println("All= " + trall +
					"  MLE= " + trmle +
					"  Max Cycles= " + MAXIT);
			int iflag2 = -1;
			// transform item parameters (BASIC 600)
			System.out.println("Reparameterized Intercepts and Slopes");
			for (int i = 0; i <= NTP - 1; i++) {
				cp[i] = cpt[0] - cpt[i+1];
				ap[i] = a[0] - a[i+1];
				if (trall == 0) {
					System.out.println("CP(" + (i + 1) + ")= " + cp[i] +
							"  AP(" + (i + 1) + ")= " + ap[i]);
				}
				else {
					// normal processing
				} // end if (trall ...
			} // end for (int i ... (BASIC 630)
			for (int nit = 1; nit <= MAXIT; nit++) {
				System.out.println("MLE Iteration= " + nit);
				// clear sums each iteration
				for (int rr = 0; rr <= NTP*2 - 1; rr++) { // rr for r
					fds[rr] = 0.0;
					pdelta[rr] = 0.0;
					for (int c = 0; c <= NTP*2 - 1; c++) {
						mtrx[rr][c] = 0.0;
					} // end for (int c ...
				} // end for (int rr ...
				for (int g = 0; g <= MAXGPS - 1; g++) {
					if (trall == 0) {
						System.out.println("Ability Group= " + (g + 1));
					}
					else {
						// normal processing
					} // end if (trall ...
					// calculate p's  (BASIC 1000)
					if (trall == 0) {
						System.out.println(
								"Evaluate Multivariate Logistic Function");
					}
					else {
						// normal processing
					} // end if (trall ...
					double d = 0;
					for (int k = 0; k <= MRC - 1; k++) {
						z[k] = cpt[k] + a[k] * theta[g];
						ez[k] = Math.exp(z[k]);
						d += ez[k];
						if (trall == 0) {
							System.out.println("CPT= " + cpt[k] +
									"  A= " + a[k] +
									"  Theta= " + theta[g]);
							System.out.println("Z(" + (k + 1) + ")= " + z[k] +
									"  EZ= " + ez[k] +
									"  D= " + d);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1040)
					for (int k = 0; k <= MRC - 1; k++) {
						p[k] = ez[k] / d;
						if (trall == 0) {
							System.out.println("P(" + (k + 1) + ")= " + p[k]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1070)
					// calculate first derivative
					if (trall == 0) {
						System.out.println("Compute First Derivatives");
					}
					else {
						// normal processing
					} // end if (trall ...
					term[0] = -r[g][1] + f[g] * p[1]; // out of for (int kk
					term[1] = -r[g][2] + f[g] * p[2];
					term[2] = theta[g] * (-r[g][1] + f[g] * p[1]);
					term[3] = theta[g] * (-r[g][2] + f[g] * p[2]);
					for (int kk = 0; kk <= NTP*2 - 1; kk++) {
						fds[kk] += term[kk];
						if (trall == 0) {
							System.out.println("Term(" + (kk + 1) +
									")= " + term[kk] +
									"  Sum First Derivative(" + (kk + 1) +
									")= " + fds[kk]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int kk ... (BASIC 1290)

					// calculate w's  (BASIC 1400)
					w[0][0] = p[1] * (1.0 - p[1]);
					w[0][1] = p[1] * p[2];
					w[1][0] = p[2] * p[1];
					w[1][1] = p[2] * (1.0 - p[2]);
					tm[0][0] = 1.0;
					tm[0][1] = theta[g];
					tm[1][0] = theta[g];
					tm[1][1] = theta[g] * theta[g];
					if (trall == 0) {
						System.out.println("Weight and Theta Matrices");
						for (int rr = 0; rr <= NTP - 1; rr++) { // rr for r
							for (int c = 0; c <= NTP - 1; c++) {
								System.out.println("R= " + (rr + 1) +
										"  C= " + (c + 1) +
										"  W= " + w[rr][c] +
										"  TM= " + tm[rr][c]);
							} // end for (int c ...
						} // end for (int rr ...
					}
					else {
						// normal processing
					} // end if (trall ... (BASIC 1490)
					// calculate second derivatives
					// Kronecker product for second derivatives (BASIC 2000)
					for (int ib = 0; ib <= NTP - 1; ib++) {
						for (int jb = 0; jb <= NTP - 1; jb++) {
							double mf = tm[ib][jb] * f[g];
							if (ib == 0 && jb == 0) {
								for (int i = 0; i <= NTP - 1; i++) {
									for (int j = 0; j <= NTP - 1; j++) {
										double temp = mf * w[i][j];
										if (i != j) {
											temp = - temp;
										}
										mtrx[i][j] += temp;
									} // end for (int j ...
								} // end for (int i ...
							}
							else {
								// normal processing
							} // end if (ib ... (BASIC 2105)
							if (ib == 0 && jb == 1) {
								for (int i = 0; i <= NTP - 1; i++) {
									for (int j = 0; j <= NTP - 1; j++) {
										double temp = mf * w[i][j];
										if (i != j) {
											temp = - temp;
										}
										mtrx[i][j+jb+1] += temp; // +1 subscript
									} // end for (int j ...
								} // end for (int i ...
							}
							else {
								// normal processing
							} // end if (ib ... (BASIC 2125)
							if (ib == 1 && jb == 0) {
								for (int i = 0; i <= NTP - 1; i++) {
									for (int j = 0; j <= NTP - 1; j++) {
										double temp = mf * w[i][j];
										if (i != j) {
											temp = - temp;
										}
										mtrx[i+ib+1][j] += temp; // +1
									} // end for (int j ...
								} // end for (int i ...
							}
							else {
								// normal processing
							} // end if (ib ... (BASIC 2145)
							if (ib == 1 && jb == 1) {
								for (int i = 0; i <= NTP - 1; i++) {
									for (int j = 0; j <= NTP - 1; j++) {
										double temp = mf * w[i][j];
										if (i != j) {
											temp = - temp;
										}
										mtrx[i+ib+1][j+jb+1] += temp; // +1 +1
									} // end for (int j ...
								} // end for (int i ...
							}
							else {
								// normal processing
							} // end if (ib ... (BASIC 2145)
						} // end for (int jb ...
					} // end for (int ib ... (BASIC 2200)
					if (trall == 0) {
						System.out.println("Information Matrix");
						for (int rr = 0; rr <= NTP*2 - 1; rr++) {
							for (int c = 0; c <= NTP*2 - 1; c++) {
								System.out.println("MTRX(" + (rr + 1) +
										"," + (c + 1) + ")= " + mtrx[rr][c]);
							} // end for (int c
						} // end for (rr ... (BASIC 2240)
					}
					else {
						// normal processing
					} // end if (trall ...
				} // end for (int g ... (BASIC 250)
				// invert matrix (BASIC 2300)
				int MBC = 3; // static final
				for (int kk = 0; kk <= NTP*2 - 1; kk++) {
					if (mtrx[0][0] <= 0.000001) {
						System.out.println("The matrix is singular, " +
								"very near singular, or indefinite");
						iflag2 = 0;
						break; // exit for (int kk ...
					}
					else {
						double sr = Math.sqrt(mtrx[0][0]); // sr for r
						for (int ik = 0; ik <= MBC - 1; ik++) {
							t[ik] = mtrx[ik+1][0] / sr;
						} // end for (ik ...
						t[NTP*2-1] = 1.0 / sr;
						for (int jk = 0; jk <= MBC - 1; jk++) {
							for (int ik = 0; ik <= MBC - 1; ik++) {
								mtrx[ik][jk] = mtrx[ik+1][jk+1] - t[ik] * t[jk];
							} // end for (ik ...
						} // end for (jk ... (BASIC 2370)
						for (int ik = 0; ik <= NTP*2 - 1; ik++) {
							mtrx[ik][NTP*2-1] = - t[ik] * t[NTP*2-1];
						} // end for (ik ...
						for (int jk = 0; jk <= MBC - 1; jk++) {
							mtrx[NTP*2-1][jk] = mtrx[jk][NTP*2-1];
						} // end for (jk ...
					} // end if (mtrx ...
				} // end for (int kk ... (BASIC 2395)
				for (int jk = 0; jk <= NTP*2 - 1; jk++) {
					for (int ik = 0; ik <= NTP*2 - 1; ik++) {
						mtrx[ik][jk] = - mtrx[ik][jk];
					} // end for (ik ...
				} // end for (jk ...
				if (trall == 0) {
					for (int jk = 0; jk <= NTP*2 - 1; jk++) {
						for (int ik = 0; ik <= NTP*2 - 1; ik++) {
							System.out.println("Row= " + (jk + 1) +
									"  Col= " + (ik + 1) +
									"  Inverse Elementmt= " + mtrx[jk][ik]);
						} // end for (ik ...
					} // end for (jk ...
				}
				else {
					// normal processing
				} // end if (trall ...
				if (iflag2 == 0) { // new line (BASIC 256)
					break; // exit for (int nit ...
				}
				else {
					// normal processing
				} // end (iflag 2 ...
				// calculate increments and new cp and ap (BASIC 2500)
				for (int rr = 0; rr <= NTP*2 - 1; rr++) {
					for (int c = 0; c <= NTP*2 - 1; c++) {
						pdelta[rr] += mtrx[rr][c] * fds[c];
					} // end for (int c ... (BASIC 2540)
					if (trall == 0 || trmle == 0) {
						System.out.println("Delta Parameter= " + pdelta[rr]);
					}
					else {
						// normal processing
					} // end if (trall ...
				} // end for (int rr ... (BASIC 2540)
				cp[0] += pdelta[0];
				cp[1] += pdelta[1];
				ap[0] += pdelta[2];
				ap[1] += pdelta[3];
				if (trall == 0 || trmle == 0) {
					System.out.println("CP(1)= " + cp[0] +
							"  AP(1)= " + ap[0]);
					System.out.println("CP(2)= " + cp[1] +
							"  AP(2)= " + ap[1]);
				}
				else {
					// normal processing
				} // end if (trall ...
				// calculate new cpt and a's (BASIC 2600)
				cpt[0] = (cp[0] + cp[1]) / 3.0;
				cpt[1] = (cp[1] - 2.0 * cp[0]) / 3.0;
				cpt[2] = (cp[0] - 2.0 * cp[1]) / 3.0;
				a[0] = (ap[0] + ap[1]) / 3.0;
				a[1] = (ap[1] - 2.0 * ap[0]) / 3.0;
				a[2] = (ap[0] - 2.0 * ap[1]) / 3.0;
				if (trall == 0 || trmle == 0) {
					for (int k = 0; k <= MRC - 1; k++) {
						System.out.println("Intercept(" + (k + 1) +
								")= " + cpt[k] +
								"Slope(" + (k + 1) +
								")= " + a[k]);
					} // end for (int k ...
				}
				else {
					// normal processing
				} // end if (trall ...
				// BASIC 290
				char converg = 'y'; // check convergence (BASIC 2800)
				// calculate change vector and obtain new estimates
				for (int rr = 0; rr <= NTP*2 - 1; rr++) {
					if (Math.abs(pdelta[rr]) > 0.05) {
						converg = 'n';
						break; // exit for (rr ...
					}
					else {
						// normal processing
					} // end if (Math.abs ...
				} // end for (int rr ... (BASIC 2820)
				if (converg == 'y') {
					System.out.println("Final Values");
					break; // exit for (nit ...
				}
				else {
					// normal processing
				} // end if (converg ...
				if (nit == MAXIT) {
					System.out.println("Reached Maximum Cycles");
				}
				else {
					// normal processing
				} // end if (nit ...
			} // end for (int nit ... (BASIC 320)
			if (iflag2 != 0) {
				// print final intercepts and slopes
				for (int k = 0; k <= MRC - 1; k++) {
					System.out.println(
							"Intercept(" + (k + 1) + ")= " + cpt[k] +
									"  Slope(" + (k + 1) + ")= " + a[k]);
				} // end for (int k ...
			}
			else {
				// normal processing
			} // end if (iflag2 ...
			System.exit(0); // last line
		} // end public static ...


/**
 * A program to estimate examinee ability parameter under
 * nominal response model for five three-category items.
 */
		static final int NITEM = 5; // Number of ITEMs

		public void nominalResponsePersonParameter() {
			double[] p = new double[MRC];
			double[][] cp = new double[NITEM][NTP];
			double[][] ap = new double[NITEM][NTP];
			double[] z = new double[MRC];
			double[] ez = new double[MRC];
			double[][] w = new double[NTP*2][NTP*2]; // matrix
			double[] u = new double[NITEM];
			double theta = 1.0; // initial value of theta
			double[][] cpt = { // read canned data
					{ 0.25, 0.50, -0.75 },
					{ 0.50, 0.25, -0.75 },
					{ -0.75, 0.50, 0.25 },
					{ -0.75, 0.50, 0.25 },
					{ 0.25, -0.75, 0.50 }
			};
			double[][] a = {
					{ -0.80, 1.00, -0.20 },
					{ 1.00, -0.80, -0.20 },
					{ 1.00, -0.20, -0.80 },
					{ 0.20, -1.00, 0.80 },
					{ -0.80, 1.00, -0.20 }
			};
			double[] irp = new double[NITEM];
			irp[0] = 3.0;
			irp[1] = 3.0;
			irp[2] = 1.0;
			irp[3] = 2.0;
			irp[4] = 2.0;
			int trall = -1; // set traces
			int trmle = -1;
			String yn = JOptionPane.showInputDialog
					("Trace all? Y/N");
			if (yn.equalsIgnoreCase("y")) {
				trall = 0;
				trmle = 0;
			}
			else{
				String ynmle = JOptionPane.showInputDialog
						("Trace MLE terms? Y/N");
				if (ynmle.equalsIgnoreCase("y")) {
					trmle = 0;
				}
				else {
					// normal processing
				} // end if (ynmle ...
			} // end if (yn ... else
			System.out.println("All= " + trall +
					"  MLE= " + trmle +
					"  Max Cycles= " + MAXIT);
			// transform item parameters (BASIC 600)
			if (trall == 0) {
				System.out.println(
						"Reparameterized Intercepts and Slopes");
			}
			else {
				// normal processing
			} // end if (trall ...
			for (int i = 0; i <= NITEM - 1; i++) {
				for (int k = 0; k <= NTP - 1; k++) {
					if (k == 0) {
						cp[i][k] = cpt[i][0] - cpt[i][1];
						ap[i][k] = a[i][0] - a[i][1];
					}
					else {
						cp[i][k] = cpt[i][0] - cpt[i][2];
						ap[i][k] = a[i][0] - a[i][2];
					} // end if (k ...
					if (trall == 0) {
						System.out.println("CP(" + (i + 1) +
								"," + (k + 1) + ")= " + cp[i][k] +
								"  AP(" + (i + 1) +
								"," + (k + 1) + ")= " + ap[i][k]);
					}
					else {
						// normal processing
					} // end if (trall ...
				} // end for (int k ...
			} // end for (int i ... (BASIC 650)
			for (int nit = 1; nit <= MAXIT; nit++) {
				System.out.println("MLE Iteration= " + nit);
				double fds = 0.0; // clear sums each iteration
				double sds = 0.0;
				for (int i = 0; i <= NITEM - 1; i++) {
					if (trall == 0) {
						System.out.println("Item= " + (i + 1));
					}
					else {
						// normal processing
					} // end if (trall ...
					// calculate p's  (BASIC 1000)
					if (trall == 0) {
						System.out.println(
								"Evaluate Multivariate Logistic Function");
					}
					else {
						// normal processing
					} // end if (trall ...
					double d = 0;
					for (int k = 0; k <= MRC - 1; k++) {
						z[k] = cpt[i][k] + a[i][k] * theta;
						ez[k] = Math.exp(z[k]);
						d += ez[k];
						if (trall == 0) {
							System.out.println("CPT= " + cpt[i][k] +
									"  A= " + a[i][k] +
									"  Theta= " + theta);
							System.out.println("Z(" + (k + 1) +
									")= " + z[k] +
									"  EZ= " + ez[k] +
									"  D= " + d);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1050)
					for (int k = 0; k <= MRC - 1; k++) {
						p[k] = ez[k] / d;
						if (trall == 0) {
							System.out.println("P(" + (k + 1) +
									")= " + p[k]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1080)
					// calculate first derivatives (BASIC 1200)
					if (trall == 0) {
						System.out.println(
								"Compute First Derivatives");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int k = 0; k <= MRC - 1; k++) {
						if (irp[i] == k  + 1.0) { // casting
							u[k] = 1.0;
						}
						else {
							u[k] = 0.0;
						} // end if (irp ...
					} // end for (int k ...
					double term = -(ap[i][0] * (u[1] - p[1]) +
							ap[i][1] * (u[2] - p[2]));
					fds += term;
					if (trall == 0) {
						System.out.println("Term= " + term +
								"  Sum First Derivative= " + fds);
					}
					else {
						// normal processing
					} // end if (trall ... (BASIC 1290)
					// calculate w's  (BASIC 1400)
					w[0][0] = p[0] * (1.0 - p[0]);
					w[0][1] = p[0] * p[1];
					w[1][0] = w[0][1];
					w[1][1] = p[1] * (1.0 - p[1]);
					w[0][2] = p[0] * p[2];
					w[1][2] = p[1] * p[2];
					w[2][2] = p[2] * (1.0 - p[2]);
					if (trall == 0) {
						System.out.println("Weights");
						System.out.println("W11= " + w[0][0] +
								"  W12=W21= " + w[0][1] +
								"  W22= " + w[1][1]);
						System.out.println("W13= " + w[0][2] +
								"  W23= " + w[1][2] +
								"  W33= " + w[2][2]);
					}
					else {
						// normal processing
					} // end if (trall ... (BASIC 1495)
					// calculate second derivatives (BASIC 2000)
					double term1 = p[0] * (ap[i][0] * ap[i][0] * w[0][1] +
							ap[i][0] * ap[i][1] * w[0][0] +
							ap[i][1] * ap[i][1] * w[0][2]) / MRC;
					double term2 = p[1] * (-ap[i][0] * ap[i][0] * w[1][1] -
							ap[i][0] * ap[i][1] * w[1][0] +
							ap[i][1] * ap[i][1] * w[1][2]) / MRC +
							ap[i][0] * ap[i][0] * w[1][1] -
							ap[i][0] * ap[i][1] * w[1][2];
					double term3 = p[2] * (ap[i][0] * ap[i][0] * w[1][2] -
							ap[i][0] * ap[i][1] * w[0][2] -
							ap[i][1] * ap[i][1] * w[2][2]) / MRC +
							ap[i][1] * ap[i][1] * w[2][2] -
							ap[i][0] * ap[i][1] * w[1][2];
					sds += term1 + term2 + term3;
					if (trall == 0) {
						System.out.println("Term1= " + term1);
						System.out.println("Term2= " + term2);
						System.out.println("Term3= " + term3 +
								"  SDS= " + sds);
					}
					else {
						// normal processing
					} // end if (trall ... (BASIC 2090)
				} // end for (int i ... (BASIC 250)
				// calculate increments and new ability (BASIC 2500)
				double tdelta = fds / sds;
				theta += tdelta;
				if (trall == 0 || trmle == 0) {
					System.out.println("Delta Parameter= " + tdelta);
					System.out.println("Ability= " + theta);
				}
				else {
					// normal processing
				} // end if (trall ...
				// check convergence (BASIC 280)
				if (Math.abs(tdelta) <= 0.05) {
					System.out.println("Reached Final Value");
					break; // exit for (nit ...
				}
				else {
					// normal processing
				} // end if (Math ...
				if (nit == MAXIT) {
					System.out.println("Reached Maximum Cycles");
				}
				else {
					// normal processing
				} // end if (nit ...
			} // end for (int nit ... (BASIC 300)
			System.out.println("Final Ability= " + theta);
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
