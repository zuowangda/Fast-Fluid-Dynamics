//#include "stdafx.h"
#include <windows.h>
#include <stdio.h>


#include "data_structure.h"
#include "utility.h"
#include "resource.h"


typedef struct {
  double number0;
  double number1;
  int command;
} ReceivedDataFormat;

typedef struct {
  int feedback;
} FFDFeedback;

typedef struct {
  double number0;
  double number1;
  int command;
}FFDSend;

typedef struct {
  int feedback;
}DymolaFeedback;

//Global Variables:
ReceivedDataFormat received;
DymolaFeedback dymolaFeedback;

//const char fFDRecvDataCh1[] = "FFDRecvDataCh1";         //Ch1: Dymola to FFD   Need define it in caption
TCHAR dymolaRecvFeedbackCh1[] = TEXT("DymolaRecvFeedbackCh1");    
//const char fFDRecvFeedbackCh2[] = "FFDRecvFeedbackCh2"; //Ch2: FFD to Dymola   Need define it in caption
TCHAR dymolaRecvDataCh2[] = TEXT("DymolaRecvDataCh2");

//******CHANNEL ONE DYMOLA TO FFD



void sendFeedbackToDymolaCh1(FFDFeedback fFDFeedback){
  HWND hSendWindow;
  HWND hRecvWindow;
  COPYDATASTRUCT sendFeedback;

  hSendWindow = GetConsoleWindow ();        //get handel of two windows
  if (hSendWindow == NULL) {
    printf("%s\n", "self handel not found");  
  }

  hRecvWindow = FindWindow(NULL, dymolaRecvFeedbackCh1);
  while(hRecvWindow == NULL){
    Sleep(500);                                     //*******set waiting time
    hRecvWindow = FindWindow(NULL, dymolaRecvFeedbackCh1);
  }
  sendFeedback.cbData = sizeof(fFDFeedback);            //set data to "package"
  sendFeedback.lpData = &fFDFeedback;
  SendMessage(hRecvWindow, WM_COPYDATA, (WPARAM)hSendWindow, (LPARAM)&sendFeedback); //send data
}

/******************************************************************************
| Receive data from Windows message
******************************************************************************/
BOOL CALLBACK receive_data_dialog(HWND hwndDlg, UINT message, 
                                  WPARAM wParam, LPARAM lParam)
{
  PCOPYDATASTRUCT pcds;

  // Hide the dialog window
  ShowWindow(hwndDlg, SW_HIDE);

  switch (message)
  {
    case WM_CLOSE:
      EndDialog(hwndDlg, 0);
      break;
    case WM_COPYDATA:
      pcds = (PCOPYDATASTRUCT)lParam;
      // If the size matches, copy the data
      if (pcds->cbData == sizeof(received)){
        memcpy_s(&received, sizeof(received), pcds->lpData, pcds->cbData);
        EndDialog(hwndDlg, 0);
      }
  }

  return FALSE;
}

/******************************************************************************
| Receive data from the other program
******************************************************************************/
void receive_data(){
  FFDFeedback ffd_feed_back;
  char msg[500];

  // Receive the data through the message
  DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_DATA), NULL, receive_data_dialog);
  
  sprintf(msg, "cosimulation.c: Received data: %d, %f, %f ", 
          received.command, received.number0, received.number1);
  ffd_log(msg, FFD_NORMAL);

  // Send feedback for successfully receive
  ffd_feed_back.feedback = 1;
  sendFeedbackToDymolaCh1(ffd_feed_back);

} // End of receive_data()

/******************************************************************************
  Read the data send from the other program
******************************************************************************/
int read_cosimulation_data(PARA_DATA *para, REAL **var)
{
  int j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL feak[1];

  // Receive data from other program
  receive_data();

  ffd_log("cosimulation.c: Sccessfully read data from the cosimulation program", FFD_NORMAL);
  return 0;
}

/******************************************************************************
|
******************************************************************************/
BOOL CALLBACK RecvDymolaFeedbackDialogCh2(HWND hwndDlg, UINT message, WPARAM wParam, LPARAM lParam) {
  PCOPYDATASTRUCT pcds; 
  ShowWindow(hwndDlg, SW_HIDE);

  switch (message){

  case WM_CLOSE:
      EndDialog(hwndDlg, 0);
      break;

    case WM_COPYDATA:
      pcds = (PCOPYDATASTRUCT)lParam;
      // If the size matches
      if (pcds->cbData == sizeof(dymolaFeedback)){
        memcpy_s(&dymolaFeedback, sizeof(dymolaFeedback), pcds->lpData, pcds->cbData);
        EndDialog(hwndDlg, 0);
      }
  }
  return FALSE;
}

int sendDataToDymolaCh2(FFDSend fFDSend){
  HWND hSendWindow;
  HWND hRecvWindow ;
  COPYDATASTRUCT sendData;

  hSendWindow = GetConsoleWindow ();        //get handel of two windows
  if (hSendWindow == NULL) {
    printf("%s\n", "self handel not found");  
    system("pause"); 
  }

  hRecvWindow = FindWindow(NULL, dymolaRecvDataCh2);
  while(hRecvWindow == NULL){
    Sleep(500);                                     //*******set waiting time
    hRecvWindow = FindWindow(NULL, dymolaRecvDataCh2);
  }
  sendData.cbData = sizeof(fFDSend);            //set data to "package"
  sendData.lpData = &fFDSend;
  SendMessage(hRecvWindow, WM_COPYDATA, (WPARAM)hSendWindow, (LPARAM)&sendData); //send data

  //run recvfeedback close dialog
  DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_FEEDBACK), NULL, RecvDymolaFeedbackDialogCh2);
  if(dymolaFeedback.feedback==1){
    return 1;
  }else {
    return 0;
  }
}


/******************************************************************************
  Write data that will be read by other program
******************************************************************************/
int write_cosimulation_data(PARA_DATA *para, REAL **var)
{

  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int check;
  FFDSend fFDSend1;

  
      //prepare data to send
  fFDSend1.number0 = 888.666;
  fFDSend1.number1 = 666.888;
  fFDSend1.command = 123456;

  check = sendDataToDymolaCh2(fFDSend1);
  if (check==0){
    ffd_log("No feedback from the other program", FFD_ERROR);
  }


  return 0;
}


