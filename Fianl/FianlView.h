
// FianlView.h: CFianlView 类的接口
//

#pragma once

typedef struct//点
{
	int x, y;
}Point;

typedef struct ENode//边表节点
{
	int xbot;//直线下端点的x坐标
	int ytop;//直线上端点的y坐标
	double fm;//斜率倒数
	struct ENode* next;//连接下一条边
}EdgeNode, * pEdgeNode;

typedef struct//边表定义
{
	pEdgeNode* base;//表头数组基地址
	int  ymin;//y=ymin为最下面一条扫描线
	int size;//表头数组长度，等于ymax+1
}EdgeTable;

typedef struct aeNode//活动边结点
{
	int ytop;//与当前扫描线yi相交的边界边的上端点的y坐标
	double xi;//当前扫描线yi与边界边的交点的x坐标
	double fm;//直线斜率倒数
	struct aeNode* next;//连接下一个边
}AENode, * pAENode;

typedef struct
{
	int x, y;
}PIXEL;

typedef struct
{
	PIXEL* pixels;//顺序表首地址
	int size;
}PixList;//用顺序表储存所有像素点


class CFianlView : public CView
{
protected: // 仅从序列化创建
	CFianlView() noexcept;
	DECLARE_DYNCREATE(CFianlView)

		// 特性
public:
	CFianlDoc* GetDocument() const;


	//直线
	bool m_isFirst = true;
	LONG m_x0 = 0;
	LONG m_y0 = 0;
	LONG m_x1 = 0;
	LONG m_y1 = 0;
	double m_r = 0;

	int m_type = 0;





	int n;//传输的数目
	int way;

	//int m_isFirst;//判断直线起点和终点
	//CPoint* point1;//存储直线的起点和终点坐标

	//CPoint* point3;// 记录圆的起点和半径
	//double m_r = 0;//圆的半径

	int x_arc[3];
	int y_arc[3];//圆弧的坐标
	//int m_isPoint;//判断圆弧的第几个点
	int m_ArcSeq = 0;//判断圆弧的第几个点

	//填充
	EdgeTable m_ET;//边表
	pAENode m_top;//活动边表的头指针
	int m_ne;//多边形边的数目
	PixList m_pixList;//像素表

	//画Bezier曲线(内含画多边形):定义法
	int m_nCtrPs;//控制点数目
	int m_nSPs;//采样点数目
	CPoint* m_pCtrPs;//控制点坐标
	CPoint* m_curve;//全部采样点坐标
	//几何作图法
	CPoint us_points[4];


	// 操作
public:
	void CFianlView::DDALine(CDC* pDC, int x0, int y0, int x1, int y1, COLORREF color);//画线段
	void CFianlView::MidPntCircle(CDC* pDC, int x0, int y0, double r, COLORREF color);//画圆
	void CFianlView::DrawArc(CDC* pDC, int x[3], int y[3], COLORREF color);//画圆弧

	int CFianlView::initAET(EdgeTable& ET, pAENode& top, int ymin);//初始化活动边表
	int CFianlView::createET(CPoint* points, int nPnts, EdgeTable& ET);//读入每条边，并插入新的节点
	int CFianlView::initET(EdgeTable& ET, int ymax, int ymin);//初始化表
	void CFianlView::intersectionPnts(int y, pAENode top, double* pXs, int& nPnts);//求当前扫描线与活动边的交点
	int CFianlView::updateAET(EdgeTable& ET, pAENode top, int yi); //剔除与当前扫描线yi不相交的边，加入与yi相交的边
	int CFianlView::polygonFill(EdgeTable& ET, pAENode top, PixList* pixList);//多边形绘制
	void CFianlView::setPixeles(int i, double start, double end, PixList* pixList);//求配对交点之间的像素点
	void CFianlView::drawPolyFilled(CDC* pDC, COLORREF color);

	void CFianlView::DrawPolygon(CDC* pDC, CPoint* m_pCtrPs, int m_nCtrPs);//绘制多边形

	//绘制Bezier曲线(定义法)
	int init();
	void CFianlView::computeCoefficients(int n, int* c);
	void CFianlView::computePoint(float t, CPoint* pt, int m_nCtrPs, CPoint* m_pCtrPs, int* c);
	void CFianlView::Bezier(CPoint* m_pCtrPs, int m_nCtrPs, int m, CPoint* m_curve);
	void CFianlView::DrawBezier(CDC* pDC, CPoint* m_curve, int m);

	//几何作图法
	void CFianlView::DrawCurve(CDC* pDC, int count);
	void CFianlView::GetPoint(CPoint* us_points);
	void CFianlView::DrawPolygon(CDC* pDC);

	// 重写
public:
	virtual void OnDraw(CDC* pDC);  // 重写以绘制该视图
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

	// 实现
public:
	virtual ~CFianlView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

	// 生成的消息映射函数
protected:
	DECLARE_MESSAGE_MAP()
public:
	//afx_msg void On_Line();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLine();
	afx_msg void OnCircle();
	afx_msg void OnArc();
	afx_msg void OnPolygon();
};

#ifndef _DEBUG  // FianlView.cpp 中的调试版本
inline CFianlDoc* CFianlView::GetDocument() const
{
	return reinterpret_cast<CFianlDoc*>(m_pDocument);
}
#endif

