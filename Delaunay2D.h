/*
 * Delaunay2D.h
 *
 *  Created on: May 18, 2013
 *      Author: nagaso
 */

#ifndef DELAUNAY2D_H_
#define DELAUNAY2D_H_

#include <cfloat>
#include <cmath>

#include <map>
#include <set>

namespace Tercel
{
  struct Vector {
    double x;
    double y;
    //double amp;//reflection signal strength

    // ======================================
    // 等価性の判定
    // ======================================
    bool operator==(const Vector& v) const
    {
      return (x == v.x && y == v.y);
    }

    // ======================================
    // 大小判定（set / mapを構築する際に使用）
    // ======================================
    bool operator<(const Vector& v) const
    {
      return x != v.x ? x < v.x : y < v.y;
    }
  };


  struct Circle
  {
    Vector center;  // 中心座標
    double radius;  // 半径
  };


  class Triangle {
  public:
    const Vector* p1, * p2, * p3;  // 頂点座標

  private:
    // 最小の頂点を返す
    // --------------------------------------
    inline const Vector* getMinVertex() const
    {
      return *p1 < *p2 && *p1 < *p3 ? p1 : *p2 < *p3 ? p2 : p3;
    }

    // まんなかの頂点を返す
    // --------------------------------------
    inline const Vector* getMidVertex() const
    {
      return *p1 < *p2 && *p2 < *p3 || *p3 < *p2 && *p2 < *p1 ? p2 :
             *p2 < *p3 && *p3 < *p1 || *p1 < *p3 && *p3 < *p2 ? p3 : p1;
    }

    // 最大の頂点を返す
    // --------------------------------------
    inline const Vector* getMaxVertex() const
    {
      return *p2 < *p1 && *p3 < *p1 ? p1 : *p3 < *p2 ? p2 : p3;
    }

  public:

    // ======================================
    // 等価性の判定
    // ======================================
    bool operator==(const Triangle& t) const
    {
      return *p1 == *t.p1 && *p2 == *t.p2 && *p3 == *t.p3 ||
             *p1 == *t.p2 && *p2 == *t.p3 && *p3 == *t.p1 ||
             *p1 == *t.p3 && *p2 == *t.p1 && *p3 == *t.p2 ||

             *p1 == *t.p3 && *p2 == *t.p2 && *p3 == *t.p1 ||
             *p1 == *t.p2 && *p2 == *t.p1 && *p3 == *t.p3 ||
             *p1 == *t.p1 && *p2 == *t.p3 && *p3 == *t.p2;
    }

    // ======================================
    // 大小判定（set / mapを構築する際に使用）
    // ======================================
    bool operator<(const Triangle& t) const
    {
      return !(*getMinVertex() == *t.getMinVertex()) ?
                 *getMinVertex() < *t.getMinVertex() :
             !(*getMidVertex() == *t.getMidVertex()) ?
                 *getMidVertex() < *t.getMidVertex() :
                 *getMaxVertex() < *t.getMaxVertex();
    }

    // ======================================
    // 他の三角形と共有点を持つか
    // ======================================
    bool hasCommonPoints(const Triangle& t) const
    {
      return *p1 == *t.p1 || *p1 == *t.p2 || *p1 == *t.p3 ||
             *p2 == *t.p1 || *p2 == *t.p2 || *p2 == *t.p3 ||
             *p3 == *t.p1 || *p3 == *t.p2 || *p3 == *t.p3;
    }
  };


  class Delaunay2d
  {
    // ======================================
    // 型定義
    // ======================================
    typedef const std::set<Vector>         ConstVertexSet;
    typedef ConstVertexSet::const_iterator ConstVertexIter;

    typedef std::set<Triangle>             TriangleSet;
    typedef std::set<Triangle>::iterator   TriangleIter;

    typedef std::map<Triangle, bool>       TriangleMap;

  private:
    // ======================================
    // 一時マップを使って重複判定
    // hashMap
    //  - Key : 三角形
    //  - Value : 重複していないかどうか
    //        - 重複していない : true
    //        - 重複している   : false
    // ======================================
    static inline void addElementToRedundanciesMap(TriangleMap* pRddcMap,
                                                   Triangle& t)
    {
      TriangleMap::iterator it = pRddcMap->find(t);
      if(it != pRddcMap->end() && it->second)
      {
        // 値を (t, true) から (t, false) に置換
        pRddcMap->erase(it);
        pRddcMap->insert(TriangleMap::value_type(t, false));
      }
      else
      {
        pRddcMap->insert(TriangleMap::value_type(t, true));
      }
    }

  public:
    static void getDelaunayTriangles(ConstVertexSet& pVertexSet,
                                     TriangleSet* triangleSet)
    {

      Triangle hugeTriangle;
      {
        // ======================================
        // 外部三角形を作る
        // ======================================
        double maxX, maxY; maxX = maxY = DBL_MIN;
        double minX, minY; minX = minY = DBL_MAX;
        for(ConstVertexIter it = pVertexSet.begin();
            it != pVertexSet.end(); ++it)
        {
          double x = it->x;
          double y = it->y;
          if(maxX < x) maxX = x; if(minX > x) minX = x;
          if(maxY < y) maxY = y; if(minY > y) minY = y;
        }

        // すべての点を包含する矩形の外接円
        double centerX  = (maxX - minX) * 0.5;      // 中心x座標
        double centerY  = (maxY - minY) * 0.5;      // 中心y座標

        double dX   = maxX - centerX;
        double dY   = maxY - centerY;
        double radius = sqrt(dX * dX + dY * dY) + 1.0;  // 半径

        Vector* p1 = new Vector;  // メモリ確保（314行目で解放）
        p1->x    = centerX - sqrt(3.0) * radius;
        p1->y    = centerY - radius;

        Vector* p2 = new Vector;  // メモリ確保（315行目で解放）
        p2->x    = centerX + sqrt(3.0) * radius;
        p2->y    = centerY - radius;

        Vector* p3 = new Vector;  // メモリ確保（316行目で解放）
        p3->x    = centerX;
        p3->y    = centerY + 2.0 * radius;

        hugeTriangle.p1 = p1;
        hugeTriangle.p2 = p2;
        hugeTriangle.p3 = p3;
      }

      triangleSet->insert(hugeTriangle);

      // --------------------------------------
      // 点を逐次添加し、反復的に三角分割を行う
      // --------------------------------------
      for(ConstVertexIter vIter = pVertexSet.begin();
          vIter != pVertexSet.end(); ++vIter)
      {
        const Vector* p = &*vIter;

        // --------------------------------------
        // 追加候補の三角形を保持する一時マップ
        // --------------------------------------
        TriangleMap rddcMap;

        // --------------------------------------
        // 現在の三角形セットから要素を一つずつ取り出して、
        // 与えられた点が各々の三角形の外接円の中に含まれるかどうか判定
        // --------------------------------------
        for(TriangleIter tIter = triangleSet->begin();
            tIter != triangleSet->end();)
        {
          // 三角形セットから三角形を取りだして…
          Triangle t = *tIter;

          // その外接円を求める。
          Circle   c;
          {
            // 三角形の各頂点座標を (x1, y1), (x2, y2), (x3, y3) とし、
            // その外接円の中心座標を (x, y) とすると、
            //     (x - x1) * (x - x1) + (y - y1) * (y - y1)
            //   = (x - x2) * (x - x2) + (y - y2) * (y - y2)
            //   = (x - x3) * (x - x3) + (y - y3) * (y - y3)
            // より、以下の式が成り立つ
            //
            // x = { (y3-y1) * (x2*x2-x1*x1+y2*y2-y1*y1)
            //     + (y1-y2) * (x3*x3-x1*x1+y3*y3-y1*y1) } / m
            //
            // y = { (x1-x3) * (x2*x2-x1*x1+y2*y2-y1*y1)
            //     + (x2-x1) * (x3*x3-x1*x1+y3*y3-y1*y1) } / m
            //
            // ただし、
            // m = 2 * {(x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)}

            double x1 = t.p1->x;  double y1 = t.p1->y;
            double x2 = t.p2->x;  double y2 = t.p2->y;
            double x3 = t.p3->x;  double y3 = t.p3->y;

            double m = 2.0*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
            double x = ((y3-y1)*(x2*x2-x1*x1+y2*y2-y1*y1)
                       +(y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;
            double y = ((x1-x3)*(x2*x2-x1*x1+y2*y2-y1*y1)
                       +(x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;

            c.center.x = x;
            c.center.y = y;

            // 外接円の半径 r は、半径から三角形の任意の頂点までの距離に等しい
            double dx   = t.p1->x - x;
            double dy   = t.p1->y - y;
            double radius = sqrt(dx * dx + dy * dy);

            c.radius = radius;
          }

          double dx = c.center.x - p->x;
          double dy = c.center.y - p->y;
          double dist = sqrt(dx * dx + dy * dy);

          // --------------------------------------
          // 追加された点が外接円内部に存在する場合、
          // その外接円を持つ三角形をリストから除外し、
          // 新たに分割し直す
          // --------------------------------------
          if(dist < c.radius)
          {
            Triangle t1;
            t1.p1 = p; t1.p2 = t.p1; t1.p3 = t.p2;
            addElementToRedundanciesMap(&rddcMap, t1);

            Triangle t2;
            t2.p1 = p; t2.p2 = t.p2; t2.p3 = t.p3;
            addElementToRedundanciesMap(&rddcMap, t2);

            Triangle t3;
            t3.p1 = p; t3.p2 = t.p3; t3.p3 = t.p1;
            addElementToRedundanciesMap(&rddcMap, t3);

            triangleSet->erase(tIter++);
          }
          else ++tIter;
        }

        // --------------------------------------
        // 一時マップのうち、重複のないものを三角形リストに追加
        // --------------------------------------
        for(TriangleMap::iterator iter = rddcMap.begin();
          iter != rddcMap.end(); ++iter)
        {
          if(iter->second) triangleSet->insert(iter->first);
        }
      }

      // --------------------------------------
      // 最後に、外部三角形の頂点を削除
      // --------------------------------------
      for(TriangleIter tIter = triangleSet->begin();
          tIter != triangleSet->end(); )
      {
        if(hugeTriangle.hasCommonPoints(*tIter))
          triangleSet->erase(tIter++);
        else ++tIter;
      }

      // 巨大三角形の頂点を解放
      delete hugeTriangle.p1;
      delete hugeTriangle.p2;
      delete hugeTriangle.p3;
    }
  };
}




#endif /* DELAUNAY2D_H_ */
