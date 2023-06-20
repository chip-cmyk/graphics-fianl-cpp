
// FianlView.cpp: CFianlView 类的实现
//
#include"cmath"
#include"iostream"
#include"string"
using namespace std;


#include "pch.h"
#include "framework.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "Fianl.h"
#endif

#include "FianlDoc.h"
#include "FianlView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#include "InputVert.h"


// CFianlView

IMPLEMENT_DYNCREATE(CFianlView, CView)

BEGIN_MESSAGE_MAP(CFianlView, CView)
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
	//ON_COMMAND(ID_2D32771, &CFianlView::On_Line)
	ON_WM_LBUTTONDOWN()
	ON_COMMAND(ID_2D32771, &CFianlView::OnLine)
	ON_COMMAND(ID_2D32772, &CFianlView::OnCircle)
	ON_COMMAND(ID_2D32773, &CFianlView::OnArc)
	ON_COMMAND(ID_2D32774, &CFianlView::OnPolygon)
	ON_COMMAND(ID_2D32775, &CFianlView::OnPolygonFill)
END_MESSAGE_MAP()

// CFianlView 构造/析构

CFianlView::CFianlView() noexcept
{
	// TODO: 在此处添加构造代码

}

CFianlView::~CFianlView()
{
}

BOOL CFianlView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CView::PreCreateWindow(cs);
}

// CFianlView 绘图

void CFianlView::OnDraw(CDC* pDC)
{
	CFianlDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	if (m_type == 1)
	{
		DDALine(pDC, m_x0, m_y0, m_x1, m_y1, RGB(255, 0, 0));
	}
	if (m_type == 2)
	{
		MidPntCircle(pDC, m_x0, m_y0, m_r, RGB(0, 0, 255));
	}
	if (m_type == 3)
	{
		DrawArc(pDC, x_arc, y_arc, RGB(0, 255, 255));
	}
	if (m_type == 4)
	{
		for (int i = 0; i < m_nCtrPs; i++)
		{
			if (i == m_nCtrPs - 1)
			{
				pDC->MoveTo(m_pCtrPs[i].x, m_pCtrPs[i].y);
				pDC->LineTo(m_pCtrPs[0].x, m_pCtrPs[0].y);

			}
			else {
				pDC->MoveTo(m_pCtrPs[i].x, m_pCtrPs[i].y);
				pDC->LineTo(m_pCtrPs[i + 1].x, m_pCtrPs[i + 1].y);
			}
			if (m_nCtrPs == n) {
				n = 0;
			}
		}
	}

	// TODO: 在此处为本机数据添加绘制代码
	//画直线
	/*switch (way)
	{
	case 1:
		DDALine(pDC, point1[0].x, point1[0].y, point1[1].x, point1[1].y, RGB(0, 0, 0));
		break;
	default:
		break;
	}*/

	//DDALine(pDC, 100, 200, 300, 400, RGB(255, 0, 0));
	//MidPntCircle(pDC, 200, 100, 50, RGB(0, 0, 255));
	//int x[3] = { 100,400,600 };
	//int y[3] = { 100,200,300 };
	//DrawArc(pDC, x, y, RGB(0, 255, 0));
	//int ymax = 0;
	//int ymin = INT_MAX;
	//m_pixList.size = 0;
	//m_pixList.pixels = new PIXEL[100000];
	//m_nCtrPs=4;
	//m_pCtrPs = new CPoint[m_nCtrPs];
	//m_pCtrPs[0] = { 300,300 };
	//m_pCtrPs[1] = { 400,300 };
	//m_pCtrPs[2] = { 400,400 };
	//m_pCtrPs[3] = { 300,400 };
	//for (int i = 0; i < m_nCtrPs; i++)
	//{
	//	if (m_pCtrPs[i].y > ymax)
	//		ymax = m_pCtrPs[i].y;
	//	if (m_pCtrPs[i].y < ymin)
	//		ymin = m_pCtrPs[i].y;
	//}
	//DrawPolygon(pDC, m_pCtrPs, m_nCtrPs);
	/*for (int i = 0; i < m_nCtrPs; i++)
	{
		if (i == m_nCtrPs - 1)
		{
			pDC->MoveTo(m_pCtrPs[i].x, m_pCtrPs[i].y);
			pDC->LineTo(m_pCtrPs[0].x, m_pCtrPs[0].y);

		}
		else {
			pDC->MoveTo(m_pCtrPs[i].x, m_pCtrPs[i].y);
			pDC->LineTo(m_pCtrPs[i + 1].x, m_pCtrPs[i + 1].y);
		}
	}*/
	//m_ne = m_nCtrPs;
	//initET(m_ET, ymax, ymin);
	//createET(m_pCtrPs, m_nCtrPs, m_ET);
	//initAET(m_ET, m_top, ymin);
	//polygonFill(m_ET, m_top, &m_pixList);
	//drawPolyFilled(pDC, RGB(120, 120, 255));
	//Invalidate();

	//m_nCtrPs = 3;
	//m_pCtrPs = new CPoint[m_nCtrPs];
	//m_pCtrPs[0] = { 100,100 };
	//m_pCtrPs[1] = { 200,200 };
	//m_pCtrPs[2] = { 300,100 };
	//m_nSPs = 500;
	//m_curve = new CPoint[m_nSPs];
	/*DrawPolygon(pDC, m_pCtrPs, m_nCtrPs);
	DrawBezier(pDC, m_curve, m_nSPs);*///定义法
	//DrawCurve(pDC, m_nCtrPs);
	//DrawPolygon(pDC);//几何作图法





}


// CFianlView 打印

BOOL CFianlView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CFianlView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加额外的打印前进行的初始化过程
}

void CFianlView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加打印后进行的清理过程
}


// CFianlView 诊断

#ifdef _DEBUG
void CFianlView::AssertValid() const
{
	CView::AssertValid();
}

void CFianlView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CFianlDoc* CFianlView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CFianlDoc)));
	return (CFianlDoc*)m_pDocument;
}
#endif //_DEBUG


// CFianlView 消息处理程序

//画直线：
void CFianlView::DDALine(CDC* pDC, int x0, int y0, int x1, int y1, COLORREF color)
{
	float dx = x1 - x0;
	float dy = y1 - y0;
	float num = abs(dx) > abs(dy) ? abs(dx) : abs(dy);
	float x = x0;
	float y = y0;
	dx = dx / (num);
	dy = dy / (num);
	for (int i = 0; i < num; i++)
	{
		pDC->SetPixel(int(x), int(y), color);
		x += dx;
		y += dy;
	}
}


//画圆：
void CFianlView::MidPntCircle(CDC* pDC, int x0, int y0, double r, COLORREF color)
{
	int x = 0;
	int y = r;
	int d = 1 - r;
	while (x <= y)
	{
		pDC->SetPixel(x0 + x, y0 + y, color);
		pDC->SetPixel(x0 + x, y0 - y, color);
		pDC->SetPixel(x0 - x, y0 + y, color);
		pDC->SetPixel(x0 - x, y0 - y, color);
		pDC->SetPixel(x0 + y, y0 + x, color);
		pDC->SetPixel(x0 + y, y0 - x, color);
		pDC->SetPixel(x0 - y, y0 + x, color);
		pDC->SetPixel(x0 - y, y0 - x, color);
		if (d < 0)
		{
			d += 2 * x + 3;
		}
		else
		{
			d += 2 * (x - y) + 5;
			y -= 1;
		}
		x += 1;
	}
}


//画圆弧：
void CFianlView::DrawArc(CDC* pDC, int x[3], int y[3], COLORREF color)
{
	double PI = 3.14;
	double dx, dy;//圆心坐标
	double a = (double)2 * (x[1] - x[0]);
	double b = (double)2 * (y[1] - y[0]);
	double c = (double)x[1] * x[1] + y[1] * y[1] - x[0] * x[0] - y[0] * y[0];
	double d = (double)2 * (x[2] - x[1]);
	double e = (double)2 * (y[2] - y[1]);
	double f = (double)x[2] * x[2] + y[2] * y[2] - x[1] * x[1] - y[1] * y[1];
	dx = (b * f - e * c) / (b * d - e * a);
	dy = (d * c - a * f) / (b * d - e * a);
	double r = sqrt((double)((x[0] - dx) * (x[0] - dx) + (y[0] - dy) * (y[0] - dy)));//半径
	double ts = atan2((double)y[0] - dy, (double)x[0] - dx);//起始弧度
	double te = atan2((double)y[2] - dy, (double)x[2] - dx);//终止弧度
	double k = (x[1] - x[0]) * (y[2] - y[1]) - (y[1] - y[0]) * (x[2] - x[1]);
	double deg(0);
	if (r < 5.08)      deg = 0.015;
	else if (r < 7.62) deg = 0.06;
	else if (r < 25.4) deg = 0.075;
	else               deg = 0.015;
	double dte = deg * 25.4 / r;
	if (k > 0)
	{
		if (x[2] > x[0] && x[1] < x[0])
		{
			te += PI * 2;
			int nCount = (int)((te - ts) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta += dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
		else if (x[2] < x[0] && x[1] < x[2])
		{
			te += PI * 2;
			int nCount = (int)((te - ts) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta += dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
		else
		{
			int nCount = (int)((te - ts) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta += dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
	}
	else if (k < 0)
	{
		if (x[2] > x[0] && x[1] < x[0])
		{
			ts += PI * 2;
			int nCount = (int)((ts - te) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta -= dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
		else if (x[2] < x[0] && x[1] < x[2])
		{
			ts += PI * 2;
			int nCount = (int)((ts - te) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta -= dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
		else
		{
			int nCount = (int)((ts - te) / dte + 0.5);
			double ta = ts;
			double xn = dx + r * cos(ts);
			double yn = dy + r * sin(ts);
			pDC->MoveTo(xn, yn);
			double ct(0), st(0);
			for (int i = 1; i <= nCount; i++)
			{
				ta -= dte;
				ct = cos(ta);
				st = sin(ta);
				xn = dx + r * ct;
				yn = dy + r * st;
				pDC->LineTo(xn, yn);
			}
		}
	}
}

//画多边形：
void CFianlView::DrawPolygon(CDC* pDC, CPoint* m_pCtrPs, int m_nCtrPs)
{
	CPen pen(PS_SOLID, 1, RGB(0, 0, 0));
	CPen* oldpen = pDC->SelectObject(&pen);
	pDC->MoveTo(m_pCtrPs[0]);
	for (int i = 1; i < m_nCtrPs; i++)
	{
		pDC->LineTo(m_pCtrPs[i]);
	}
	pDC->LineTo(m_pCtrPs[0]);
	pDC->SelectObject(oldpen);
}


//多边形填充：
int CFianlView::initAET(EdgeTable& ET, pAENode& top, int ymin)
{
	pEdgeNode p = ET.base[ymin];
	pAENode q, rear;
	top = new AENode;
	if (!top) return 0;
	top->next = NULL;
	rear = top;
	while (p)
	{
		q = new AENode;
		rear->next = q;
		q->fm = p->fm;
		q->xi = p->xbot;
		q->ytop = p->ytop;
		q->next = NULL;
		rear = rear->next;
		p = p->next;
	}
	return 1;
}

int CFianlView::initET(EdgeTable& ET, int ymax, int ymin)
{
	ET.base = new pEdgeNode[ymax + 1];
	if (!ET.base)
	{
		printf("存储空间分配失败！");
		return 0;
	}
	for (int i = 0; i < ymax + 1; i++)
	{
		ET.base[i] = new EdgeNode;
		ET.base[i]->xbot = 0;
		ET.base[i]->ytop = 0;
		ET.base[i]->next = NULL;
		ET.base[i]->fm = 1.0;
	}
	ET.size = ymax + 1;
	ET.ymin = ymin;
	return 1;
}

int CFianlView::createET(CPoint* points, int nPnts, EdgeTable& ET)
{
	int i, j, b, t, deltay, xbot, ybot, ytop;
	double fm;
	pEdgeNode p, q, r;
	for (i = 0; i < nPnts; i++)
	{
		j = (i + 1) % nPnts;
		b = (points[i].y < points[j].y) ? i : j;
		t = (i == b) ? j : i;
		xbot = points[b].x;
		ybot = points[b].y;
		ytop = points[t].y;
		deltay = ytop - ybot;
		fm = (points[t].x - xbot) / (1.0 * deltay);
		q = new EdgeNode;
		q->ytop = ytop;
		q->xbot = xbot;
		q->fm = fm;
		q->next = NULL;
		p = ET.base[ybot];
		while (p->next != NULL)
		{
			if (p->next->xbot > q->xbot)
			{
				q->next = p->next;
				p->next = q;
				break;
			}
			p = p->next;
		}
		if (q->next == NULL)
		{
			p->next = q;
		}
	}
	return 1;
}

void CFianlView::intersectionPnts(int y, pAENode top, double* pXs, int& nPnts)
{
	int i;
	bool lastTag, side;
	double x;
	nPnts = 0;
	pAENode p = top->next;
	while (p)
	{
		x = p->xi;
		if (p->ytop == y) side = false;
		else side = true;
		pXs[nPnts] = x;
		if (nPnts > 0 && pXs[nPnts] == pXs[nPnts - 1] && lastTag != side)
		{
			nPnts--;
		}
		else
		{
			nPnts++;
		}
		lastTag = side;
		p = p->next;
	}
}

int CFianlView::updateAET(EdgeTable& ET, pAENode top, int yi)
{
	pAENode p, q, r;
	pEdgeNode pEN;
	q = top;
	while (q)
	{
		if (q->next != NULL && (q->next->ytop) <= yi)
		{
			q->next = q->next->next;
		}
		q = q->next;
	}
	p = top->next;
	while (p)
	{
		p->xi += p->fm;
		p = p->next;
	}
	pEN = ET.base[yi];
	r = top;
	while (pEN)
	{
		if (pEN->ytop > yi)
		{
			q = new AENode;
			q->fm = pEN->fm;
			q->xi = pEN->xbot;
			q->ytop = pEN->ytop;
			q->next = NULL;
			while (r->next != NULL)
			{
				if ((r->next->xi) > (q->ytop))
				{
					q->next = r->next;
					r->next = q;
					break;
				}
				r = r->next;
			}
			if (q->next == NULL)
			{
				r->next = q;
			}
		}
		pEN = pEN->next;
	}
	return 1;
}

int CFianlView::polygonFill(EdgeTable& ET, pAENode top, PixList* pixList)
{
	int i = ET.ymin, n = ET.size;
	int j, k;
	int nPnts;
	double* pXs = new double[m_ne];
	while (i < n)
	{
		intersectionPnts(i, m_top, pXs, nPnts);
		j = 0;
		k = 1;
		while (k < nPnts)
		{
			setPixeles(i, pXs[j], pXs[k], pixList);
			j = k + 1;
			k = j + 1;
		}
		i++;
		if (i < n)
		{
			updateAET(m_ET, m_top, i);
		}
	}
	delete pXs;
	return 1;
}

void CFianlView::setPixeles(int i, double start, double end, PixList* pixList)
{
	int n, xstart, xend;
	double t;
	if (start > end)
	{
		t = start;
		start = end;
		end = t;
	}
	xstart = (int)(start + 0.5);
	xend = (int)(end + 0.5);
	while (xstart <= xend)
	{
		PIXEL* pixel = new PIXEL;
		pixel->x = xstart;
		pixel->y = i;
		pixList->pixels[pixList->size] = *pixel;
		pixList->size++;
		xstart += 1;
	}
}

void CFianlView::drawPolyFilled(CDC* pDC, COLORREF color)
{
	for (int i = 0; i < m_pixList.size + 1; i++)
	{
		int x = m_pixList.pixels[i].x;
		int y = m_pixList.pixels[i].y;
		pDC->SetPixel(x, y, color);
	}
}

//绘制Bezier曲线
int CFianlView::init()
{
	//将空间的首地址赋给m_pCtrPs和m_curve指针（注意在适当时候释放相应的存储空间）
	m_nCtrPs = 0;
	m_nSPs = n * 100;
	m_pCtrPs = new CPoint[n];  //最多10个点
	m_curve = new CPoint[m_nSPs];   //假定400个
	if (!m_pCtrPs)
	{
		return 0;
	}
	//	CPoint *m_curve; // （全部采样点坐标）
	return 1;
}
void CFianlView::computeCoefficients(int n, int* c)//计算第k个点的系数c[k]
{//n 为控制点数目，c 为存储空间的首地址，存储内容为系数
	int k, i;
	for (k = 0; k <= n; k++)
	{
		c[k] = 1;
		for (i = n; i >= k + 1; i--) /*求 c[k]=n*(n-1)…(k+1) */
			c[k] *= i;
		for (i = n - k; i >= 2; i--) /*求 c[k]/(n-k)!*/
			c[k] /= i;
	}
}

void CFianlView::computePoint(float t, CPoint* pt, int m_nCtrPs, CPoint* m_pCtrPs, int* c)
{
	//计算 Bezier 曲线上参数为 t 的点
	//pt 为所求点，nCtrPs 为控制点数目，pCtrPs 为存储控制点坐标的空间首地址
	int i;
	float n = m_nCtrPs - 1;
	float blend;
	float t1 = 1 - t; //基函数的值
	pt->x = 0.0;
	pt->y = 0.0;
	for (i = 0; i <= n; i++)
	{
		blend = c[i] * powf(t, i) * powf(t1, n - i);
		pt->x += m_pCtrPs[i].x * blend; //求 x(t)
		pt->y += m_pCtrPs[i].y * blend; //求 y(t)
	}
}

void CFianlView::Bezier(CPoint* m_pCtrPs, int m_nCtrPs, int m, CPoint* m_curve)
{
	//m 个采样点，结果保存在 curve 所指的数组里面
	int i;
	int* pC = (int*)malloc(m_nCtrPs * sizeof(int)); //分配系数的存储空间
	computeCoefficients(m_nCtrPs - 1, pC);
	for (i = 0; i <= m; i++)
		computePoint(i / (float)m, &m_curve[i], m_nCtrPs, m_pCtrPs, pC);
	free(pC);
}

void CFianlView::DrawBezier(CDC* pDC, CPoint* m_curve, int m)//绘制Bezier曲线
{
	Bezier(m_pCtrPs, m_nCtrPs, m, m_curve);
	CPen pen(PS_SOLID, 2, RGB(255, 0, 0));
	CPen* oldpen = pDC->SelectObject(&pen);
	pDC->MoveTo(m_curve[0]);
	for (int i = 1; i < m; i++)
	{
		pDC->LineTo(m_curve[i]);
	}
	//pDC->SelectObject(oldpen);
}

//Bezier几何作图法
void CFianlView::DrawCurve(CDC* pDC, int count)
{
	CPoint points[20][20];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			points[i][j] = { 0,0 };
		}
	}
	double t = 0.0;
	for (int i = 0, j = 0; i < count; i++)
	{
		points[i][j] = m_pCtrPs[i];
	}
	DrawPolygon(pDC, m_pCtrPs, m_nCtrPs);
	for (t = 0.0; t < 1.0; t += 0.0001)
	{
		for (int j = 1; j < count; j++)
		{
			for (int i = 0; i < count - j; i++)
			{
				if (points[i][j - 1].x > points[i + 1][j - 1].x) {   //该曲线收到points[i][j-1]的影响为t
					points[i][j].x = points[i][j - 1].x - abs(points[i][j - 1].x - points[i + 1][j - 1].x) * t;
				}
				else {
					points[i][j].x = points[i][j - 1].x + abs(points[i][j - 1].x - points[i + 1][j - 1].x) * t;
				}

				if (points[i][j - 1].y > points[i + 1][j - 1].y) {
					points[i][j].y = points[i][j - 1].y - abs(points[i][j - 1].y - points[i + 1][j - 1].y) * t;
				}
				else {
					points[i][j].y = points[i][j - 1].y + abs(points[i][j - 1].y - points[i + 1][j - 1].y) * t;
				}
			}

		}
		pDC->SetPixel((int)points[0][count - 1].x, (int)points[0][count - 1].y, RGB(0, 0, 255));
	}
}

void CFianlView::GetPoint(CPoint* us_points)
{
	for (int i = 0; i < m_nCtrPs; i++)
	{
		us_points[i] = m_pCtrPs[i];
	}
}


void CFianlView::DrawPolygon(CDC* pDC)
{
	CPen pen;
	CPen* oldpen;
	pen.CreatePen(PS_SOLID, 2, RGB(200, 200, 200));
	oldpen = pDC->SelectObject(&pen);
	GetPoint(m_pCtrPs);
	pDC->MoveTo(us_points[0].x, us_points[0].y);
	pDC->Ellipse(us_points[0].x - 5, us_points[0].y - 5, us_points[0].x + 5, us_points[0].y + 5);
	for (int i = 1; i < 4; i++)
	{
		pDC->LineTo(us_points[i].x, us_points[i].y);
		pDC->Ellipse(us_points[i].x - 5, us_points[i].y - 5, us_points[i].x + 5, us_points[i].y + 5);
	}
	pDC->SelectObject(oldpen);
	pen.DeleteObject();
}

//void CFianlView::On_Line()
//{
//	// TODO: 在此添加命令处理程序代码
//	point1 = new CPoint[2];
//	m_isFirst = 0;
//	way = 1;
//	Invalidate();
//}



void CFianlView::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	if (m_type == 1 || m_type == 2)
	{
		if (m_isFirst)// m_isFirst 为 CxxView 的成员，第一次按下鼠标左键
		{
			m_isFirst = false;
			m_x0 = point.x;
			m_y0 = point.y;
		}
		else
		{
			m_isFirst = true;
			m_x1 = point.x;
			m_y1 = point.y;
			int dx = abs(m_x1 - m_x0);
			int dy = abs(m_y1 - m_y0);
			m_r = sqrt(dx * dx + dy * dy);
			Invalidate(true);//刷新屏幕
		}
	}
	if (m_type == 3)
	{
		if (m_ArcSeq == 0)
		{
			x_arc[0] = point.x;
			y_arc[0] = point.y;
			//AfxMessageBox(_T("圆弧圆弧圆弧"));
			m_ArcSeq++;
		}
		else if (m_ArcSeq == 1)
		{
			x_arc[1] = point.x;
			y_arc[1] = point.y;
			m_ArcSeq++;

		}
		else if (m_ArcSeq == 2)
		{
			x_arc[2] = point.x;
			y_arc[2] = point.y;
			m_ArcSeq = 0;
			Invalidate(true);//刷新屏幕
		}
	}
	if (m_type == 4)
	{
		if (m_nCtrPs < n)
		{
			m_pCtrPs[m_nCtrPs] = point;
			m_nCtrPs++;  //采样点增
		}
		Invalidate(true);
	}
	CView::OnLButtonDown(nFlags, point);
}


void CFianlView::OnLine()
{
	// TODO: 在此添加命令处理程序代码
	m_type = 1;

}


void CFianlView::OnCircle()
{
	// TODO: 在此添加命令处理程序代码
	m_type = 2;

}


void CFianlView::OnArc()
{
	// TODO: 在此添加命令处理程序代码
	m_type = 3;
}


void CFianlView::OnPolygon()
{
	// TODO: 在此添加命令处理程序代码
	m_type = 4;
	InputVert VertDlg;
	if (IDOK == VertDlg.DoModal())
	{
		n = VertDlg.m_InputVert;
	}
	init();

}


void CFianlView::OnPolygonFill()
{
	// TODO: 在此添加命令处理程序代码
	m_type = 5;
}
