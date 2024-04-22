#ifndef VPT_QT_COMMON_H
#define VPT_QT_COMMON_H

#include <QColor> 
#include <QIcon>

#include <osg/Matrix> 
#include <osg/Vec3>
#include <osg/Vec4>

namespace vpt{

class ColorIcon : public QIcon{
protected: 
public: 
	ColorIcon(const QColor& color, const QRect& rect){
		QPixmap pmap(rect.width(),rect.height()); 
		pmap.fill(color);
		addPixmap(pmap);
		//paint(painter, rect); 
	}
}; 

inline QColor osg2QtColor(const osg::Vec4& color){
	//return QColor(std::min(255,(int)(color.x()*255)), std::min(255,(int)(color.y()*255)), std::min(255,(int)(color.z()*255))); 
	return QColor((int)(color.x()*255),
		(int)(color.y()*255), (int)(color.z()*255)); 

}

inline osg::Vec4 Qt2osgColor(const QColor& color){
	return osg::Vec4(color.redF(),color.greenF(),color.blueF(),color.alphaF()); 
}

}
#endif
