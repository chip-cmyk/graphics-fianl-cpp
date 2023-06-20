#pragma once
#include "afxdialogex.h"


// InputVert 对话框

class InputVert : public CDialogEx
{
	DECLARE_DYNAMIC(InputVert)

public:
	InputVert(CWnd* pParent = nullptr);   // 标准构造函数
	virtual ~InputVert();

	// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG1 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	int m_InputVert;
};
