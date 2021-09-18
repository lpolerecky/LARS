# LARS

<img src="man/figures/LARS-logo.png">

**Look@Rates** - a Matlab script for calculating *substrate assimilation rates* in cells from their isotopic composition determined by nanoSIMS.


WORK IN PROGRESS.

## Input data in a spreadsheet

Input data needs to be organized in a spreadsheet. A template spreadsheet is available in the folder *templates*.

The required data include the parameters listed in the table below. Each row corresponds to an individual cell. Although the explanation assumes carbon isotopes (i.e., <sup>13</sup>C atom fraction), the values can be entered for isotopes of any other element (e.g., <sup>15</sup>N, <sup>18</sup>O, etc.).

| row  | explanation |
|-----:|:------------|
| **t**    | Duration of the SIP incubation. A value in *h* will yield rates in *h*<sup>-1</sup>.|
| **x**    | Best estimate of the <sup>13</sup>C atom fraction of the measured cell, *x*(<sup>13</sup>C). The value is determined from the nanoSIMS measurement. It can be exported by analyzing the nanoSIMS data with Look@NanoSIMS.|
| **dx** | Error (uncertainty) of the <sup>13</sup>C atom fraction of the measured cell, Δ*x*(<sup>13</sup>C). Value determined from the nanoSIMS measurement. It can be exported by analyzing the nanoSIMS data with Look@NanoSIMS. |
| **xSeff** | <sup>13</sup>C atom fraction of the effective carbon source, *x*(<sup>13</sup>C)<sub>S,eff</sub> (see Eq. 6). Ideally, the value is obtained from a direct measurement. Alternatively, the value can be calculated from the known concentration of substrate in the medium and the known concentration and amount of added labelled substrate. |
| **x_ini** | Initial <sup>13</sup>C atom fraction of the cell, *x*(<sup>13</sup>C)<sub>i</sub>. Ideally, the value is determined from the measurement of control cells. If not available, the value corresponding to the natural abundance can be used (e.g., *x*(<sup>13</sup>C)<sub>i</sub> = 0.011, *x*(<sup>15</sup>N)<sub>i</sub> = 0.0037, *x*(<sup>18</sup>O)<sub>i</sub> = 0.002, etc.) |
| **rho** | Carbon density of the measured cell, &#961;. A value reflecting the average carbon content per &#956;m<sup>3</sup> of the cell (in fmol C &#956;m<sup>-3</sup>).  |
| **avgVcell** | Average biovolume of the cell from the same species as the measured cell, &#10216;V&#10217; (in &#956;m<sup>3</sup>). The values of &#961; and &#10216;V&#10217; are used to calculate the average C content of the cell as &#10216;C&#10217; = &#961;&#183;&#10216;V&#10217;. Note that the calculated value must represent the C content of the cell averaged across the entire cell cycle. Thus, it can be determined by measuring bulk C content and cell counts for a population of cells with perfectly unsynchronized cell cycles, but not from similar measurements for a population with partially synchronized cells. Alternatively, it can be obtained by measuring the biovolume or carbon content of dividing cells (V<sub>max</sub> or C<sub>max</sub>, respectively), as explained in the manuscrit (Approach 3A).|
| **Vcell** | Biovolume of the measured cell, V (in &#956;m<sup>3</sup>). This value is relevant if Approach B or C is used to calculate the rate (see explanation for **dVcell**). Specifically, the carbon content of the measured cell is calculated as C = &#961;&#183;V. Ideally, V is calculated from the size and shape of the cell determined from the same nanoSIMS image as the atom fraction, x. Note that V must be such that the cell cycle stage, calculated from the biovolume as s = V/(V<sub>max</sub>/2)-1, ends up between 0 and 1. In this calculation, V<sub>max</sub> is calculated from the average biovolume as V<sub>max</sub> = &#10216;V&#10217; &#183; 2 &#183; ln(2). If the entered combination of values of &#10216;V&#10217; and V yields s<0 or s>1, the calculation will automatically change V to the value of &#10216;V&#10217;. |
| **dVcell** | Error (uncertainty) of the measured cell’s biovolume, &#916;V (in &#956;m<sup>3</sup>). |

WORK IN PROGRESS.