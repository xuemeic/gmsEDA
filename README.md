# Electrodermal Activity Decomposition
`gmsEDA.m` is the main function. For usage, see `example1.m` or `example2.m`

## Data
`eda_5m32s_23ptcp_shift.mat` contains EDA signal of 23 subjects. 
```
load("eda_5m32s_23ptcp_shift.mat")
eda = eda_5m32s_23ptcp_shift.eda;
```
`eda` is a matrix with 23 columns. Each column is the EDA signal of a subject of 5 min and 32 seconds long, measured at 4Hz. **Video starts at 6th seconds**.

### video content
- 0:00 – 0:28:  Neutral video, scene of people walking on a busy city street
- 0:29 – 0:59:  Neutral video, scene of two people working together on a laptop
- 1:00 – 2:14:  Black screen with text ``Please complete distractor task now''
- 2:15 – 2:42:  Positive video, baby attempting but failing to drink water from a hose
	- 2:22: first laugh moment
	- 2:35: Baby smiles
- 2:43 – 3:14:  Positive video, cat staring at camera, wiggling tongue
	- Video is consistent throughout (no unique moments)
- 3:15 – 4:28:  Black screen with text ``Please complete distractor task now''
- 4:29 – 4:57:  Negative video, skateboarder falling and breaking arm
	- 4:34:  Moment of fall (doesn’t look overly disturbing)
	- 4:39:  Broken arm clearly shown (very disturbing) 
- 4:58 – 5:25:  &Negative video, animal trainer has arm chomped by alligator
	- 5:17:  Alligator chomps on arm and begins to roll with arm in mouth
