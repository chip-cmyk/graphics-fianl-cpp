// InputVert.cpp: 实现文件
//

#include "pch.h"
#include "Fianl.h"
#include "afxdialogex.h"
#include "InputVert.h"


// InputVert 对话框

IMPLEMENT_DYNAMIC(InputVert, CDialogEx)

InputVert::InputVert(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG1, pParent)
	, m_InputVert(0)
{

}

InputVert::~InputVert()
{
}

void InputVert::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_InputVert);
}


BEGIN_MESSAGE_MAP(InputVert, CDialogEx)
END_MESSAGE_MAP()


// InputVert 消息处理程序
