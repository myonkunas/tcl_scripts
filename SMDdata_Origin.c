/* The Origin C program used to process the SMD data. The input data contains 
 * Following information: Time, z, force,(z,force)... from 4 independent runs. 
 * Processing process includes:
 * (1). Convert Time information to position of reference point
 * (2). Integrate the force into work
 * (3). Average the position of ion during 10 ps windows
 * (4). Calculate the average work and squred average work in the window
 * (5). Reconstruct PMF.
 * Adapted from Automation.c example file from Origin
 * by Zhanwu Liu, University of Pittsburgh. Sep 02, 2003
 */

////////////////////////////////////////////////////////////////////
//
#include <origin.h>
//
////////////////////////////////////////////////////////////////////

void SMD_Processing (string strFileName)
{
	//The length of data, SMD data output every 20 steps, total 150,000 steps
	int iXlen = 7500;
	//Extract just the name of the file passed, from the full file name with path
	int iLen=strFileName.GetLength();
	int iSlash=strFileName.ReverseFind('\\');
	int iDot = strFileName.ReverseFind('.');
	string strSampleName = strFileName.Mid(iSlash+1,iDot-iSlash-1);
	printf("Processing file : %s ..........", strSampleName);

	//Declare root folder of PE for the current project
	Folder fldRootFolder = Project.RootFolder;

	//Create subfolder with same name as file
	Folder fldSubFolder = fldRootFolder.AddSubfolder(strSampleName);

	//Make this subfolder active
	fldSubFolder.Activate();

	//Get file path to the current project
	string strProjectPath;
	strProjectPath = Project.GetPath();

	//Creat a worksheet using custom template
	//Generate a template worksheet by myself
	Worksheet wksData;
	//string strWksTemplate = strProjectPath + "SMDtemplate.otw";
	string strWksTemplate = "z:\SMDtemplate.otw";
	int nOptionW = CREATE_VISIBLE_SAME; //enum of different status
	bool bRetW = wksData.Create(strWksTemplate,nOptionW);

	//Declare datasets in worksheet to copy data from file
	Dataset dsT(wksData,0);
	Dataset dsZ1(wksData,1);
	Dataset dsF1(wksData,2);
	Dataset dsZ2(wksData,3);
	Dataset dsF2(wksData,4);
	Dataset dsZ3(wksData,5);
	Dataset dsF3(wksData,6);
	Dataset dsZ4(wksData,7);
	Dataset dsF4(wksData,8);

	//open file and read data into worksheet
	stdioFile ffDataFile;
	ffDataFile.Open(strFileName, file::modeRead);

	//Read first two lines as headers.
	//First line contains: velocity(A/ps), beta(1/kBT), starting point, direction (1 or -1);
	//Second line contains: T Z1 F1 Z2 F2 Z3 F3 Z4 F4
	
        string strParameters, strHeaderLine, strData;
	ffDataFile.ReadString(strParameters);
	ffDataFile.ReadString(strHeaderLine);

	string str1 = strParameters.GetToken(0);
	string str2 = strParameters.GetToken(1);
	string str3 = strParameters.GetToken(2);
	string str4 = strParameters.GetToken(3);

	//Convert string to numbers 
	double velocity = atof(str1);
	double beta = atof(str2);
	double startPoint = atof(str3);
	int direction = atoi(str4);


	//set dataset length appropriately, defined in the first part of program
	dsT.SetSize(iXlen); 
        dsZ1.SetSize(iXlen);
        dsF1.SetSize(iXlen);
        dsZ2.SetSize(iXlen);
        dsF2.SetSize(iXlen);
 	dsZ3.SetSize(iXlen);
	dsF3.SetSize(iXlen);      
        dsZ4.SetSize(iXlen);
        dsF4.SetSize(iXlen);


	//Loop thru and read all data points, read in all the data. 
	for (int ii = 0; ii < iXlen; ii++) {
		ffDataFile.ReadString(strData);
		string str0 = strData.GetToken(0);
           	string str1 = strData.GetToken(1);
		string str2 = strData.GetToken(2);
		string str3 = strData.GetToken(3);
		string str4 = strData.GetToken(4);
           	string str5 = strData.GetToken(5);
		string str6 = strData.GetToken(6);
		string str7 = strData.GetToken(7);
		string str8 = strData.GetToken(8);

		dsT[ii] =  atof(str0);
                dsZ1[ii] = atof(str1);
                dsF1[ii] = atof(str2);
	        dsZ2[ii] = atof(str3);
                dsF2[ii] = atof(str4);
                dsZ3[ii] = atof(str5);
                dsF3[ii] = atof(str6);
                dsZ4[ii] = atof(str7);
                dsF4[ii] = atof(str8);
		printf("processing line %f \n",dsT[ii]);
	}
}

