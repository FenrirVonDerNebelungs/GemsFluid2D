// GemsFluid2DWin.cpp : Defines the entry point for the application.
//


#include "GemsFluid2DWin.h"
#include "GenImage.h"
#include "Fluid2D.h"
#include "Test/Test.h"

#define MAX_LOADSTRING 100
#define G_ID_SIZE 1001
#define G_ID_SIM_FREE_RUN 1002
#define G_ID_CUDA_LOOPS 1003
#define G_ID_VISC 1004
#define G_ID_JAC_LOOPS 1005
#define G_ID_JAC_F_LOOPS 1006
#define G_ID_RUN_BUTTON 103

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
HWND g_fluid_hWnd;
HWND g_setup_L1_hWnd;/*size*/
HWND g_setup_1_hWnd;
HWND g_setup_L2_hWnd;/*sim free run loops*/
HWND g_setup_2_hWnd;
HWND g_setup_L3_hWnd;/*internal cuda cycles*/
HWND g_setup_3_hWnd;
HWND g_setup_L4_hWnd;/*viscosity exponent*/
HWND g_setup_4_hWnd;
HWND g_setup_L5_hWnd;/*jacobi loops*/
HWND g_setup_5_hWnd;
HWND g_setup_L6_hWnd; /*jacobi force loops*/
HWND g_setup_6_hWnd;
HWND g_setup_button_hWnd;
GenImage* g_pGenImage = nullptr; //class used to draw images
Fluid2D* g_pFluid2D = nullptr; //class used to run the fluid simulation
Test* g_pTest = nullptr; //class used to run the fluid simulation
bool g_state_run = false;

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
VOID                Release();
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_GEMSFLUID2DWIN, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_GEMSFLUID2DWIN));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {    
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        if (g_state_run) {
            if (g_pFluid2D->runSim()) {
                InvalidateRect(g_fluid_hWnd, nullptr, TRUE);
                g_pGenImage = g_pFluid2D->getCurrentImage();
            }
        }
    }

    Release();
    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_GEMSFLUID2DWIN));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_GEMSFLUID2DWIN);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
	/*g_pFluid2D = new Fluid2D();
	s_WH grid_wh = g_pFluid2D->getGridWidthHeight();*/
	//g_pFluid2D->applyForce();
	//g_pFluid2D->initFileOutput();
	//g_pFluid2D->launchCUDA();
	//g_pFluid2D->releaseFileOutput();


	//g_pTest = new Test();
	//grid_wh = g_pTest->getGridWidthHeight();
    //g_pTest->runTest();
    //g_pGenImage = g_pTest->getTestImage();//new GenImage(grid_wh.width, grid_wh.height);
    //g_pGenImage = g_pFluid2D->getCurrentImage();
   hInst = hInstance; // Store instance handle in our global variable
   //s_WH blown_grid_wh = g_pTest->getBlownWidthHeight();
   int window_width = 512+100;//grid_wh.width + 100;//blown_grid_wh.width + 6;
   int window_height = 512+100;//grid_wh.height + 100;//blown_grid_wh.height + 6;
   g_fluid_hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, window_width, window_height, nullptr, nullptr, hInstance, nullptr);//CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!g_fluid_hWnd)
   {
      return FALSE;
   }
   g_setup_L1_hWnd = CreateWindowW(L"Static", L"Size 1 to 5 ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 40, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_1_hWnd=CreateWindowW(L"Edit", L"4", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 40, 200, 25, g_fluid_hWnd, (HMENU)G_ID_SIZE, nullptr, nullptr);
   g_setup_L2_hWnd = CreateWindowW(L"Static", L"Sim Free Run ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 80, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_2_hWnd = CreateWindowW(L"Edit", L"9", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 80, 200, 25, g_fluid_hWnd, (HMENU)G_ID_SIM_FREE_RUN, nullptr, nullptr);
   g_setup_L3_hWnd = CreateWindowW(L"Static", L"CUDA loops ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 120, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_3_hWnd = CreateWindowW(L"Edit", L"12", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 120, 200, 25, g_fluid_hWnd, (HMENU)G_ID_CUDA_LOOPS, nullptr, nullptr);
   g_setup_L4_hWnd = CreateWindowW(L"Static", L"Viscosity e- ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 160, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_4_hWnd = CreateWindowW(L"Edit", L"5", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 160, 200, 25, g_fluid_hWnd, (HMENU)G_ID_VISC, nullptr, nullptr);
   g_setup_L5_hWnd = CreateWindowW(L"Static", L"Jacobi Loops ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 200, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_5_hWnd = CreateWindowW(L"Edit", L"24", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 200, 200, 25, g_fluid_hWnd, (HMENU)G_ID_JAC_LOOPS, nullptr, nullptr);
   g_setup_L6_hWnd = CreateWindowW(L"Static", L"Jacobi F Loops ", WS_VISIBLE | WS_CHILD | SS_LEFT, 20, 240, 200, 25, g_fluid_hWnd, nullptr, nullptr, nullptr);
   g_setup_6_hWnd = CreateWindowW(L"Edit", L"64", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_NUMBER, 210, 240, 200, 25, g_fluid_hWnd, (HMENU)G_ID_JAC_F_LOOPS, nullptr, nullptr);
   g_setup_button_hWnd = CreateWindowW(L"Button", L"Run", WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON, 40, 380, 230, 30, g_fluid_hWnd, (HMENU)G_ID_RUN_BUTTON, nullptr, nullptr);
   if (!g_setup_L1_hWnd || !g_setup_1_hWnd || !g_setup_button_hWnd)
       return FALSE;

   ShowWindow(g_fluid_hWnd, nCmdShow);
   UpdateWindow(g_fluid_hWnd);

   return TRUE;
}
VOID Release()
{
    /*if (g_pTest != nullptr)
    {
        delete g_pTest;
        g_pTest = nullptr;
    }*/

    if(g_pFluid2D != nullptr)
    {
        g_pFluid2D->releaseFileOutput();
        delete g_pFluid2D;
        g_pFluid2D = nullptr;
	}
}
//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE: Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//  WM_LBUTTONDOWN: - mouse button down, stores the location of the mouse click
//  WM_LBUTTONUP: - mouse button up, applies the force based on the mouse movement, or restarts simulation if mouse drag was too small
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    static int last_mouse_button_down_x = -1;
	static int last_mouse_button_down_y = -1;

    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Parse the menu selections:
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            case G_ID_RUN_BUTTON:
            {
                wchar_t buffer1[256], buffer2[256], buffer3[256], buffer4[256], buffer5[256], buffer6[256];
                GetWindowTextW(g_setup_1_hWnd, buffer1, 256);
                GetWindowTextW(g_setup_2_hWnd, buffer2, 256);
                GetWindowTextW(g_setup_3_hWnd, buffer3, 256);
                GetWindowTextW(g_setup_4_hWnd, buffer4, 256);
                GetWindowTextW(g_setup_5_hWnd, buffer5, 256);
                GetWindowTextW(g_setup_6_hWnd, buffer6, 256);
                int val1 = _wtoi(buffer1);
                int sim_frames = _wtoi(buffer2);
                int cuda_loops = _wtoi(buffer3);
                int visc_exp = _wtoi(buffer4);
                int jac_loops = _wtoi(buffer5);
                int jac_force_loops = _wtoi(buffer6);
                int block_size = 1;
                int pow_size_cnt = 0;
                do {
                    block_size *= 2;
                    pow_size_cnt++;
                } while (block_size < 32 && pow_size_cnt < val1);
                double viscosity_nu = 1.0e-5;
                if (visc_exp <= 3)
                    viscosity_nu = 1.0e-3;
                else if (visc_exp == 4)
                    viscosity_nu = 1.0e-4;
                else if (visc_exp >= 6)
                    viscosity_nu = 1.0e-6;
                g_pFluid2D = new Fluid2D(
                    60,/*mouse max delta*/
                    5, /*mouse min delta*/
                    sim_frames,
                    cuda_loops,
                    6,/*max force frame duration */
                    3,/*force decay frames*/
                    0.1,/*force decay factor */
                    300.0,/*max allowed force*/
                    3,/*max dye frames duration*/
                    1.0e-3,/*delta_t*/
                    1.0e-3,/*delta_x*/
                    viscosity_nu,
                    block_size,/* num blocks */
                    16, /* num threads */
                    2, /* jacobi min blocks side dim*/
                    4, /* jacobi min threads side dim*/
                    jac_loops, /* number of jacobi loops per plate in stack */
                    jac_force_loops); /* number of jacobi loops '           ' when force is present */
                g_state_run = true;
                if (g_setup_1_hWnd != nullptr 
                    && g_setup_L1_hWnd != nullptr
                    && g_setup_L2_hWnd!= nullptr
                    && g_setup_2_hWnd!=nullptr
                    && g_setup_L3_hWnd != nullptr
                    && g_setup_3_hWnd != nullptr
                    && g_setup_L4_hWnd != nullptr
                    && g_setup_4_hWnd != nullptr
                    && g_setup_L5_hWnd != nullptr
                    && g_setup_5_hWnd != nullptr
                    && g_setup_L6_hWnd != nullptr
                    && g_setup_6_hWnd != nullptr
                    && g_setup_button_hWnd != nullptr) {
                    DestroyWindow(g_setup_L1_hWnd);
                    DestroyWindow(g_setup_1_hWnd);
                    DestroyWindow(g_setup_L2_hWnd);
                    DestroyWindow(g_setup_2_hWnd);
                    DestroyWindow(g_setup_L3_hWnd);
                    DestroyWindow(g_setup_3_hWnd);
                    DestroyWindow(g_setup_L4_hWnd);
                    DestroyWindow(g_setup_4_hWnd);
                    DestroyWindow(g_setup_L5_hWnd);
                    DestroyWindow(g_setup_5_hWnd);
                    DestroyWindow(g_setup_L6_hWnd);
                    DestroyWindow(g_setup_6_hWnd);
                    DestroyWindow(g_setup_button_hWnd);
                    break;
                }
            }
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }                
        }
        break;
    case WM_PAINT:
        {
        if (g_pGenImage != nullptr && g_state_run) {
            PAINTSTRUCT ps;
            BITMAPINFO* bmi = g_pGenImage->getBitmapInfo();
            s_WH wh = g_pGenImage->getWidthHeight();
            unsigned char* pImageData = g_pGenImage->getImageData();

            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            if (pImageData != nullptr) {
                SetDIBitsToDevice(hdc,
                    0, 0,               // xDest, yDest
                    wh.width, wh.height,     // nWidth, nHeight
                    0, 0,             // xSrc, ySrc
                    0,                // Start scan line
                    wh.height,          // number of scan lines
                    pImageData, // pointer to the array of color data
                    bmi,
                    DIB_RGB_COLORS);
            }
            EndPaint(hWnd, &ps);
        }

        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_LBUTTONDOWN: {
        last_mouse_button_down_x = LOWORD(lParam);
        last_mouse_button_down_y = HIWORD(lParam);
        //g_pGenImage = g_pFluid2D->handleMouse();//g_pTest->handleMouse();
        //if(g_pTest->getMessage()!=nullptr)
        //    MessageBox(hWnd, g_pTest->getMessage(), L"Active", MB_OK);
		//InvalidateRect(hWnd, nullptr, TRUE);/*invalidates the entire client area and causes a WM_PAINT message to be sent to the window procedure
                                              //NULL means entire client area, true to erase background*/
        //MessageBox(hWnd, L"Left mouse button clicked", L"Mouse Click", MB_OK);
        break;
    }
    case WM_LBUTTONUP: {
		int x_release = LOWORD(lParam);
        int y_release = HIWORD(lParam);
        if(g_state_run)
           g_pFluid2D->handleMouseSweep(x_release, y_release, last_mouse_button_down_x, last_mouse_button_down_y);
		last_mouse_button_down_x = -1;
		last_mouse_button_down_y = -1;
        //MessageBox(hWnd, L"Left mouse button released", L"Mouse Click", MB_OK);
		break;
    }
    case WM_KEYDOWN: {
        switch (wParam)
        {
            case VK_SPACE:
                if(g_state_run)
                    g_pFluid2D->advanceSim();
                break;
            case VK_ESCAPE:
                if(g_state_run)
                    g_pFluid2D->resetSim();
                break;
            case VK_OEM_PLUS:
                if(g_state_run)
                    g_pFluid2D->toggleAddDyeWithForce();
                break;
            case VK_OEM_MINUS:
                if (g_state_run)
                    g_pFluid2D->toggleAddDyeWithoutForce();
        }
        break;
    }

    case WM_CHAR: {
        switch (wParam)
        {
        case 'p':
            if(g_state_run){
                g_pFluid2D->setDisplayP();
                g_pGenImage = g_pFluid2D->getCurrentImage();
                InvalidateRect(hWnd, nullptr, TRUE);
                break;
            }
            case 'd': 
            if(g_state_run){
                g_pFluid2D->setDisplayDye();
                g_pGenImage = g_pFluid2D->getCurrentImage();
                InvalidateRect(hWnd, nullptr, TRUE);
                break;
            }
            case 's':
            if(g_state_run){
                 g_pFluid2D->toggleSlidingWallActive();
                if (g_pFluid2D->slidingWallActive())
                    MessageBox(hWnd, L"activating sliding wall", L"State", MB_OK);
                else
                    MessageBox(hWnd, L"turning off sliding wall", L"State", MB_OK);
                break;
            }
            case 'w':
            if(g_state_run){
                if (g_pFluid2D->initFileOutput())
                    MessageBox(hWnd, L"Dumping output to disk", L"File Out", MB_OK);
                else
                    MessageBox(hWnd, L"init dump failed, file output may already be active", L"File Out", MB_OK);
                break;
            }
            default:
                break;
        }
            break;
    }
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
