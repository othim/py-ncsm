The format of the interaction files is:

-------------------------------------------------------------------------------
<meta_data>
n   l   np  lp  s   j   mt  ME
DATA:
0	0	0	0	0	0	-1	-7.651247774229011
0	0	1	0	0	0	-1	-3.234632983474619
0	0	2	0	0	0	-1	1.472983482380104
0	0	3	0	0	0	-1	5.851416257753321
0	0	4	0	0	0	-1	9.453923188728981
0	0	5	0	0	0	-1	11.80509728323482
...
-------------------------------------------------------------------------------

- The interaction matric element is defined as 
  
  ME = <n l s j t mt | V_12 | np lp s j t mt> (MeV)

  where n,l,s,j,t,mt are the quantum numbers for the two nucleon system in 
  relative coordinates. V_12 is the two-nucleon potential.

- Conventions:

  - The potential matrix elements follow the Machleidt convention.

  - The radial harmonic oscillator wave functions R_{nl}(p) have the (-1)^n 
    factor in momentum space.
