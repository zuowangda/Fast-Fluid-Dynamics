//#include "stdafx.h"
#include <windows.h>
#include <stdio.h>


#include "data_structure.h"
#include "utility.h"
#include "resource.h"
#include "cosimulation.h"

//Global Variables:
ReceivedData received;
ReceivedCommand command_received;



//const char fFDRecvDataCh1[] = "FFDRecvDataCh1";         //Ch1: Dymola to FFD   Need define it in caption
TCHAR command_window_other[] = TEXT("DymolaRecvFeedbackCh1");    
//const char fFDRecvFeedbackCh2[] = "FFDRecvFeedbackCh2"; //Ch2: FFD to Dymola   Need define it in caption
TCHAR data_window_other[] = TEXT("DymolaRecvDataCh2");



/******************************************************************************
| Receive the command from the other program
******************************************************************************/
BOOL CALLBACK receive_command(HWND hwndDlg, UINT message, WPARAM wParam, 
                              LPARAM lParam) 
{
  PCOPYDATASTRUCT pcds; 

  // Hide message window
  ShowWindow(hwndDlg, SW_HIDE);

  switch(message)
  {
    case WM_CLOSE:
      EndDialog(hwndDlg, 0);
      break;
    case WM_COPYDATA:
      pcds = (PCOPYDATASTRUCT)lParam;
      // If the size matches, copy the data
      if (pcds->cbData == sizeof(command_received))
      {
        memcpy_s(&command_received, sizeof(command_received), pcds->lpData, pcds->cbData);
        EndDialog(hwndDlg, 0);
      }
  }
  return FALSE;

} // End of receive_command()

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

  ffd_log("cosimulation.c: Sccessfully read data from the cosimulation program", 
           FFD_NORMAL);
  return 0;
}

/******************************************************************************
| Receive data from the other program
******************************************************************************/
void receive_data()
{
  SentCommand command_sent;
  char msg[500];

  // Receive the data through the message
  DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_DATA), NULL, receive_data_dialog);
  
  sprintf(msg, "cosimulation.c: Received data: %d, %f, %f ", 
          received.command, received.number0, received.number1);
  ffd_log(msg, FFD_NORMAL);

  // Send feedback for successfully receive
  command_sent.feedback = 1;
  send_command(command_sent);

} // End of receive_data()



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
} //End of receive_data_dialog()



/******************************************************************************
| Send commands to the other programs
******************************************************************************/
void send_command(SentCommand command_sent){
  HWND hSendWindow;
  HWND hRecvWindow;
  COPYDATASTRUCT data_package;

  // Get handel of two windows
  hSendWindow = GetConsoleWindow ();        
  if (hSendWindow == NULL)
    ffd_log("cosimulation.c: self handel not found", FFD_ERROR);

  // Check if the other program is ready (created relative console window) 
  // for receiving the message
  hRecvWindow = FindWindow(NULL, command_window_other);

  // Continue to check the status if the other program is ready
  while(hRecvWindow == NULL){
    // Wait
    Sleep(500);
    // Check it again
    hRecvWindow = FindWindow(NULL, command_window_other);
  }

  // Set the command to a package
  data_package.cbData = sizeof(command_sent);
  data_package.lpData = &command_sent;

  // Send the command
  SendMessage(hRecvWindow, WM_COPYDATA, (WPARAM)hSendWindow, (LPARAM)&data_package);
}

/******************************************************************************
| Send data to the other program
******************************************************************************/
int send_data(SentData data_sent){
  HWND hSendWindow;
  HWND hRecvWindow ;
  COPYDATASTRUCT data_package;

  //get handel of two windows
  hSendWindow = GetConsoleWindow ();        
  if (hSendWindow == NULL) {
    ffd_log("cosimulation.c: Self handel not found", FFD_ERROR);
    system("pause"); 
  }

  // Get the window handle of the other program
  hRecvWindow = FindWindow(NULL, data_window_other);
  while(hRecvWindow == NULL){
    // Wait
    Sleep(500);
    // Check again
    hRecvWindow = FindWindow(NULL, data_window_other);
  }

  // Set data to a package
  data_package.cbData = sizeof(data_sent);
  data_package.lpData = &data_sent;

  // Send the package
  SendMessage(hRecvWindow, WM_COPYDATA, (WPARAM)hSendWindow, (LPARAM)&data_package);

  // Check if the data package have been received
  DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_FEEDBACK), NULL, receive_command);

  if(command_received.feedback==1)
    return 0;
  else
    return 1;
} // End of send_data()


/******************************************************************************
  Write data that will be read by other program
******************************************************************************/
int write_cosimulation_data(PARA_DATA *para, REAL **var)
{

  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  SentData data_sent;
  
  //prepare data to send
  data_sent.number0 = 888.666;
  data_sent.number1 = 666.888;
  data_sent.command = 11111;

  if (!send_data(data_sent))
    ffd_log("cosimulation.c: Successfully sent data to the other program", FFD_NORMAL);
  else
    ffd_log("cosimulation.c: No feedback from the other program", FFD_ERROR);

  return 0;
}


