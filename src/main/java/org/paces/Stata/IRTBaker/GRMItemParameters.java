package org.paces.Stata.IRTBaker;

import javax.swing.*;

/**
 * Created by billy on 9/27/15.
 */
public class GRMItemParameters {

/**
 * Maximum likelihood estimation of item parameters for a graded
 * response item with four ordered response categories.
 */
		static final int MRC = 4; // nuMber of Response Categories
		static final int MBC = 3; // nuMber of Broundary Curves
		static final int MAXIT = 4; // MAXimum number of ITerations
		static final int MAXGPS = 7; // number of theta Points on Scale
		public GRMItemParameters(String[] args) {
			this.gradedResponseItemParameters();
			this.gradedResponsePersonParameter();
		}

		public void gradedResponseItemParameters() {
			double[] pstr = new double[MRC+1]; // for BASIC 1200
			double[] w = new double[MRC+1]; // for BASIC 1200
			double[] pk = new double[MRC];
			double[] pobs = new double[MRC];
			double[] sfd = new double[MRC];
			double[] t = new double[MRC];
			double[] pdelta = new double[MRC];
			double[] b = new double[MBC];
			double[] bk = new double[MRC];
			double[][] mtrx = new double[MRC][MRC];
			double[][] tk = new double[MRC][MRC];
			// read canned data
			double a = 1.0; // initial value of common slope
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
					{ 691.0, 242.0, 61.0, 6.0 },
					{ 500.0, 341.0, 136.0, 23.0 },
					{ 309.0, 382.0, 242.0, 67.0 },
					{ 159.0, 341.0, 341.0, 159.0 },
					{ 67.0, 242.0, 382.0, 309.0 },
					{ 23.0, 136.0, 341.0, 500.0 },
					{ 6.0, 61.0, 242.0, 691.0 },
			};
			double[] cpt = new double[MBC];
			cpt[0] = 1.8;
			cpt[1] = 0.0;
			cpt[2] = -1.8;
			int trall = -1; // establish trace of program
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
			int iflag = -1;
			for (int nit = 1; nit <= MAXIT; nit++) {
				// clear sums, vectors, and matrix
				for (int rr = 0; rr <= MRC - 1; rr++) { // rr for r
					sfd[rr] = 0.0;
					for (int c = 0; c <= MRC - 1; c++) {
						mtrx[rr][c] = 0.0;
					} // end for (int c ...
				} // end for (int rr ...
				System.out.println("MLE Iteration= " + nit);
				for (int g = 0; g <= MAXGPS - 1; g++) {
					System.out.println("Ability Level= " + (g + 1));
					// calc p star
					pstr[0] = 1.0;
					w[0] = 0.0;
					pstr[MRC] = 0.0;
					w[MRC] = 0.0;
					if (trall == 0) {
						System.out.println("Boundard Curve Pstar's");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int y = 0; y <= MBC - 1; y++) {
						double dev = cpt[y] + a * theta[g];
						pstr[y+1] = 1.0 / (1.0 + Math.exp(-dev));
						w[y+1] = pstr[y+1] * (1.0 - pstr[y+1]);
						if (trall == 0) {
							System.out.println("Intercept(" + (y + 1) + ")= " + cpt[y] +
									"  A= " + a +
									"  Theta(" + (g + 1) + ")= " + theta[g]);
							System.out.println("Dev= " + dev +
									"  Pstar(" + (y + 1) + ")= " + pstr[y+1] +
									"  W(" + (y + 1) + ")= " + w[y+1]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int y ... (BASIC 1280)
					// calc p observed
					if (trall == 0) {
						System.out.println("Pk and P observed");
					}
					else {
						// normal processing
					}
					for (int k = 0; k <= MRC - 1; k++) {
						if (k == 0) {
							pk[k] = 1.0 - pstr[k+1];
						}
						else if (k == MRC - 1) {
							pk[MRC-1] = pstr[MRC-1];
						}
						else {
							pk[k] = pstr[k] - pstr[k+1];
						} // end if (k ...
						pobs[k] = r[g][k] / f[g];
						if (trall == 0) {
							System.out.println("K= " + (k + 1) +
									"  PK= " + pk[k] +
									"  PObs= " + pobs[k]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1360)
					// calculate first derivative terms for intercepts
					if (trall == 0) {
						System.out.println("Intercept Terms");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int y = 0; y <= MBC - 1; y++) {
						double term1 = - pobs[y] / pk[y] + pobs[y+1] / pk[y+1];
						double term2 = f[g] * w[y+1] * term1;
						sfd[y] += term2;
						if (trall == 0) {
							System.out.println("Y= " + (y + 1) +
									"  Term1= " + term1 +
									"  Term2= " + term2 +
									"  First Derivative= " + sfd[y]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int y ... (BASIC 1460)
					// calculate first derivative terms for slope
					double sumls = 0.0;
					if (trall == 0) {
						System.out.println("Common Slope Terms");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int k = 0; k <= MRC - 1; k++) {
						double term1 = w[k] - w[k+1];
						double term2 = (pobs[k] / pk[k]) * term1;
						sumls += term2;
						if (trall == 0) {
							System.out.println("K= " + (k + 1) +
									"  Term1= " + term1 +
									"  Term2= " + term2 +
									"  SumLS= " + sumls);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int k ... (BASIC 1560)
					double term3 = f[g] * theta[g] * sumls;
					sfd[MRC-1] = sfd[MRC-1] + term3;
					if (trall == 0) {
						System.out.println("Term3= " + term3 +
								"  First Derivative Term= " + sfd[MRC-1]);
					}
					else {
						// normal processing
					} // end if (trall ...
					// calculate intercept terms for information matrix
					if (trall == 0) {
						System.out.println(
								"Intercept Terms in Information Matrix");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int rr = 0; rr <= MBC - 1; rr++) { // rr for r
						for (int c = 0; c <= MBC - 1; c++) {
							if (rr == c) { // diagonal intercept term
								double term1 = 1 / pk[c] + 1 / pk[c+1];
								double term2 = f[g] * w[c+1] * w[c+1] * term1;
								mtrx[rr][c] += term2;
								if (trall == 0) {
									System.out.println("R= " + (rr + 1) +
											"  C= " + (c + 1) +
											"  Term1= " + term1 +
											"  Term2= " + term2);
									System.out.println(
											"Matrix Element= " + mtrx[rr][c]);
								}
								else {
									// normal processing
								} // end if (trall ...
							}
							else {
								double term1 = f[g] * w[rr+1] * w[c+1];
								double term2 = - term1 * (1.0 / pk[c]); // switched
								if (rr > c) {
									term2 = - term1 * (1.0 / pk[rr]);
								}
								else {
									// normal processing
								} // end if (rr ...
								mtrx[rr][c] += term2;
								if (Math.abs(rr - c) > 1) {
									mtrx[rr][c] = 0.0;
								}
								else {
									// normal processing
								} // end if (Math.abs ...
								if (trall == 0) {
									System.out.println("R= " + (rr + 1) +
											"  C= " + (c + 1) +
											"  Term1= " + term1 +
											"  Term2= " + term2);
									System.out.println(
											"Matrix Element= " + mtrx[rr][c]);
								}
								else {
									// normal processing
								} // end if (trall ...
							} // end if (rr ...
						} // end for (int c ... (BASIC 1730)
					} // end for (int rr ... (BASIC 1730)
					// calculate slope terms for information matrix
					if (trall == 0) {
						System.out.println(
								"Slope Terms in Information Matrix");
					}
					else {
						// normal processing
					} // end if (trall ...
					int rr = MRC - 1; // rr for r
					for (int c = 0; c <= MRC - 1; c++) {
						if (trall == 0) {
							System.out.println("R= " + (rr + 1) +
									"  C= " + (c + 1));
						}
						else {
							// normal processing
						} // end if (trall ...
						double term1 = 0.0; // initialize
						double term2 = 0.0;
						if (c != MRC - 1) {
							term1 = (w[c] - w[c+1]) / pk[c] -
									(w[c+1] - w[c+2]) / pk[c+1];
							term2 = - f[g] * w[c+1] * theta[g] * term1;
						}
						else {
							term2 = 0.0;
							for (int k = 0; k <= MRC -1; k++) {
								term1 = (w[k] - w[k+1]) * (w[k] - w[k+1]) *
										1.0 / pk[k];
								term2 += term1;
								if (trall == 0) {
									System.out.println("K= " + (k + 1) +
											"  Term1= " + term1 +
											"  Term2= " + term2);
								}
								else {
									// normal processing
								} // end if (trall ...
							} // end for (int k ...
							term2 = term2 * f[g] * theta[g] * theta[g];
						} // end if (c ...
						mtrx[rr][c] += term2;
						if (trall == 0) {
							System.out.println("Term1= " + term1 +
									"  Term2= " + term2 +
									"  Matrix Element= " + mtrx[rr][c]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (int c ... (BASIC 1930)
					// put slope row derivatives in last column
					for (int c = 0; c <= MRC - 1; c++) {
						mtrx[c][rr] = mtrx[rr][c];
					} // end for (int c ...
					if (trall == 0) {
						System.out.println("Information Matrix");
					}
					else {
						// normal processing
					} // end if (trall ...
					for (int row = 0; row <= MRC - 1; row++) { // row for r
						for (int c = 0; c <= MRC - 1; c++) {
							if (trall == 0) {
								System.out.println("R= " + (row + 1) +
										"  C= " + (c + 1) +
										"  Matrix Element= " + mtrx[row][c]);
							}
							else {
								// normal processing
							} // end if (trall ...
						} // end for (int c ...
					} // end for (int r ...
				} // end for (int g ... (BASIC 270)
				// invert information matrix
				for (int kk = 0; kk <= MRC - 1; kk++) {
					if (mtrx[0][0] <= 0.000001) {
						System.out.println("The matrix is singular, " +
								"very near singular, or indefinite");
						iflag = 0;
						break; // exit for (int kk ...
					}
					else {
						double sr = Math.sqrt(mtrx[0][0]); // sr for r
						for (int ik = 0; ik <= MBC - 1; ik++) {
							t[ik] = mtrx[ik+1][0] / sr;
						} // end for (ik ...
						t[MRC-1] = 1.0 / sr;
						for (int jk = 0; jk <= MBC - 1; jk++) {
							for (int ik = 0; ik <= MBC - 1; ik++) {
								mtrx[ik][jk] = mtrx[ik+1][jk+1] - t[ik] * t[jk];
							} // end for (ik ...
						} // end for (jk ... (BASIC 2170)
						for (int ik = 0; ik <= MRC - 1; ik++) {
							mtrx[ik][MRC-1] = - t[ik] * t[MRC-1];
						} // end for (ik ...
						for (int jk = 0; jk <= MBC - 1; jk++) {
							mtrx[MRC-1][jk] = mtrx[jk][MRC-1];
						} // end for (jk ...
					} // end if (mtrx ...
				} // end for (int kk ... (BASIC 2195)
				for (int jk = 0; jk <= MRC - 1; jk++) {
					for (int ik = 0; ik <= MRC - 1; ik++) {
						mtrx[ik][jk] = - mtrx[ik][jk];
					} // end for (ik ...
				} // end for (jk ...
				if (trall == 0) {
					for (int jk = 0; jk <= MRC - 1; jk++) {
						for (int ik = 0; ik <= MRC - 1; ik++) {
							System.out.println("Row= " + (jk + 1) +
									"  Col= " + (ik + 1) +
									"  Inversion Elementmt= " + mtrx[jk][ik]);
						} // end for (ik ...
					} // end for (jk ...
				}
				else {
					// normal processing
				} // end if (trall ...
				// BASIC 290
				char converg = 'y';
				if (iflag == 0) {
					break; // exit for (int nit ...
				}
				else {
					// calculate change vector and obtain new estimates
					for (int jk = 0; jk <= MRC - 1; jk++) {
						pdelta[jk] = 0.0;
					} // end for (jk ...
					for (int jk = 0; jk <= MRC - 1; jk++) {
						for (int ik = 0; ik <= MRC - 1; ik++) {
							pdelta[jk] += mtrx[jk][ik] * sfd[ik];
						} // end for (ik ...
						if (trall == 0 || trmle == 0) {
							System.out.println("Delta Parameter(" + (jk + 1) +
									")= " + pdelta[jk]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (jk ...
					// may need big delta protection here
					// obtain new item parameter vlaues
					for (int ik = 0; ik <= MBC - 1; ik++) {
						cpt[ik] += pdelta[ik];
						if (trall == 0 || trmle == 0) {
							System.out.println("Intercept(" + (ik + 1) +
									")= " + cpt[ik]);
						}
						else {
							// normal processing
						} // end if (trall ...
					} // end for (ik ...
					a += pdelta[MRC-1];
					if (trall == 0 || trmle == 0) {
						System.out.println("Common Slope= " + a);
					}
					else {
						// normal processing
					} // end if (trall ...
					// calculate convergence criterion
					for (int k = 0; k <= MRC - 1; k++) {
						if (Math.abs(pdelta[k]) > 0.05) {
							converg = 'n';
							break; // exit for (k ...
						}
						else {
							// normal processing
						} // end if (Math.abs ...
					} // end for (k ...
				} // end if (iflag ...
				if (converg == 'y') {
					System.out.println("Converg=y"); // new line
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
			} // end for (int nit ... (BASIC 330)
			if (iflag != 0) {
				// print item parameter estimates
				for (int y = 0; y <= MBC - 1; y++) {
					b[y] = - cpt[y] / a;
					System.out.println("Intercept(" + (y + 1) +
							")= " + cpt[y] +
							"Boundary Curve B(" + (y + 1) +
							")= " + b[y]);
				} // end for (int y ...
				System.out.println("Common Slope= " + a);
				for (int k = 0; k <= MRC - 1; k++) {
					if (k == 0) {
						bk[k] = b[0];
					}
					else if (k == MRC - 1) {
						bk[k] = b[MBC - 1];
					}
					else{
						bk[k] = (b[k] + b[k-1]) / 2.0;
					}
					System.out.println("Response Category(" + (k + 1) +
							") Difficulty= " + bk[k]);
				} // end for (k ...
				System.out.println("Item Discrimination= " + a);
			}
			else {
				// normal processing
			} // end if (iflag ...
			System.exit(0); // last line
		} // end public static ...

/**
 * A program to estimate ability under graded response model for
 * five four-response test items.
 */
		static final int NITEM = 5; // Nunber of ITEMs
		public void gradedResponsePersonParameter() {
			double[] a = new double[NITEM]; // read canned data
			a[0] = 0.85;
			a[1] = 1.00;
			a[2] = 0.90;
			a[3] = 1.10;
			a[4] = 1.40;
			double[][] cpt = {
					{ 1.72, 0.00, -1.72 },
					{ 1.20, 0.20, -1.50 },
					{ 2.10, 1.10, -0.70 },
					{ -0.60, -0.80, -1.20 },
					{ 1.40, 0.70, -1.70 }
			};
			double[] u = { 3.0, 2.0, 3.0, 4.0, 3.0 };
			double theta = 0.4; // initial value of theta
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
						("Trace MLE? Y/N");
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
			double fds = 0.0; // first derivative sum
			double sds = 0.0; // second derivative sum
			System.out.println("Initial Theta= " + theta);
			for (int nit = 1; nit <= MAXIT; nit++) {
				System.out.println("MLE Iteration= " + nit);
				for (int i = 0; i <= NITEM - 1; i++) {
					double term4 = 0.0;
					for (int k = 0; k <= MRC - 1; k++) {
						double k1 = k + 1.0; // casting
						if (u[i] == k1) {
							char ptsflg = 'n';
							// calc p's and w's
							int km1 = k - 1;
							double pstrkm1 = 1.0;
							double wkm1 = 0.0;
							double devkm1 = 0.0; // initialize
							if (km1 != -1){
								devkm1 = cpt[i][km1] + a[i] * theta;
								pstrkm1 = 1.0 / (1.0 + Math.exp(-devkm1));
								wkm1 = pstrkm1 * (1.0 - pstrkm1);
							}
							else {
								// normal processing
							} // end if (km1 ...
							double pstrk = 0.0;
							double wk = 0.0;
							double dev = 0.0; // initialize
							if (k != MRC -1) {
								dev = cpt[i][k] + a[i] * theta;
								pstrk = 1.0 / (1.0 + Math.exp(-dev));
								wk = pstrk * (1.0 - pstrk);
							}
							else {
								// normal processing
							} // end if (k ...
							double pk = pstrkm1 - pstrk;
							if (trall == 0) {
								System.out.println("Item(" + (i + 1) + ")" +
										"  U(I)= " + u[i]);
								System.out.println("K-1=" + k +
										"  Dev(K-1)= " + devkm1 +
										"  Pstar(K-1)= " + pstrkm1 +
										"  W(K-1)= " + wkm1);
								System.out.println("K=" + (k + 1) +
										"  Dev(K)= " + dev +
										"  Pstar(K)= " + pstrk +
										"  W(K)= " + wk);
							}
							else {
								// normal processing
							} // end if (trall ...
							if (pk < 0.0001 || pk > 0.9999) {
								ptsflg = 'y';
								break; // exit for (int k ...
							}
							else {
								// normal processing
							} // end if (pk ...
							// calc first and second derivative terms
							// first derivative
							double fterm = (wkm1 - wk) / pk;
							fds += fterm;
							if (trall == 0) {
								System.out.println("Fterm= " + fterm +
										"  First Derivative Sum= " + fds);
							}
							else {
								// normal processing
							} // end if (trall ...
							// second derivative
							double term1 = -wk * (1.0 - 2.0 * pstrk) / pk;
							double term2 = wkm1 * (1.0 - 2.0 * pstrkm1) / pk;
							double term3 = - (wkm1 - wk) * (wkm1 - wk) /
									(pk * pk);
							term4 += term1 + term2 + term3;
							sds += a[i] * a[i] * term4;
							if (trall == 0) {
								System.out.println("Term1= " + term1 +
										"  Term2= " + term2 +
										"  Term3= " + term3);
								System.out.println("Term4= " + term4 +
										"  Second Derivative Sum= " + sds);
							}
							else {
								// normal processing
							} // end if (trall ...
						}
						else {
							// normal processing
						} // end if (u[i] ...
					} // end for (int k ... (BASIC 230)
				} // end for (int i ... (BASIC 240)
				double delta = fds / sds;
				theta -= delta;
				if (trall == 0 || trmle == 0) {
					System.out.println("Theta= " + theta +
							"  Delta= " + delta);
				}
				else {
					// normal processing
				} // end if (trall ...
				if (Math.abs(delta) < 0.05) {
					break; // exit for (nit ...
				}
				else {
					// normal processing
				} // end if (Math.abs ...
				if (nit == MAXIT) {
					System.out.println("Maximum Cycles Reached");
				}
				else {
					// normal processing
				} // end if (nit ...
			} // end for (int nit ... (BASIC 280)
			System.out.println("Theta Estimate= " + theta);
			System.exit(0); // last line
		} // end public static ...
	} // end public class ...
