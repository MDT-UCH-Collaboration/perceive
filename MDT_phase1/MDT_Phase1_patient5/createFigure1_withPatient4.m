% Create Figure 1 for patient 4

jsonFiles = 'pt4_0819_r_survey.json'; 
js = jsondecode(fileread(jsonFiles));

[outTABLE] = extractSurveyData(js,'r',1);

%%
dbdir = 'C:\Users\Admin\Documents\Github\perceive\MDT_phase1\MDT_Phase1_patient5'
mdt_p1_p5figure(dbdir,4)