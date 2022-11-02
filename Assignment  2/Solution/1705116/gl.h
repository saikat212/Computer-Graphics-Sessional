#ifndef GL_H
#define GL_H


#include "bitmap_image.hpp"

#include<bits/stdc++.h>
using namespace std;

class NullBuffer : public std::streambuf
{
public:
	int overflow(int c) { return c; }
} nullBuffer;

ostream nullStream(&nullBuffer);

template<typename T>
class Gl
{
	std::stack<Mat4<T>> m_stack;
	ostream * stage1;
	ostream * stage2;
	ostream * stage3;
	ostream * z_out;
	Mat4<T> m_view,m_proj;
	vector<Triangle<T>> m_triangles;

	public:
	Gl():stage1( & nullStream),
			 stage2( & nullStream),
			 stage3( &nullStream),
			 z_out( &nullStream)
	{
		m_stack.push(Mat4<T>::identity());
	}
	void setStage1(ostream &s)
	{
		stage1 = &s;
	}
	
	void setStage2(ostream &s)
	{
		stage2 = &s;
	}
	
	void setStage3(ostream &s)
	{
		stage3 = &s;
	}

	void setZOut(ostream &s)
	{
		z_out = &s;
	}

	void lookAt(
							Vec3<T> eye,
							Vec3<T> look,
							Vec3<T> up
						)
	{
		Vec3<T> l = look-eye;
		l.normalize();
		Vec3<T> r = l.cross(up);
		r.normalize();
		Vec3<T> u = r.cross(l);
		Mat4<T> tran = Mat4<T>::translate(-eye);
		
		Mat4<T> rot = Mat4<T>(); // zero
		for(int i=0;i<3;i++)	rot[0][i] = r[i];
		for(int i=0;i<3;i++)	rot[1][i] = u[i];
		for(int i=0;i<3;i++)	rot[2][i] = -l[i];
		rot[3][3] = T(1);

		m_view = rot * tran;
		
	}

	void perspective(T fovY, T aspectRatio, T near, T far)
	{
		T fovX = fovY * aspectRatio;
		T t = near * tan(fovY*PI/T(360));
		T r = near * tan(fovX*PI/T(360));

		m_proj = Mat4<T>(); // zero matrix
		m_proj[0][0] = near/r;
		m_proj[1][1] = near/t;
		m_proj[2][2] = -(far+near)/(far-near);
		m_proj[2][3] = -T(2)*far*near/(far-near);
		m_proj[3][2] = -T(1);
	}

	void scale(Vec3<T>  &v)
	{
		m_stack.top() = m_stack.top() * Mat4<T>::scale(v);
	}
	void translate(Vec3<T> const &v)
	{
		m_stack.top() = m_stack.top() * Mat4<T>::translate(v);
	}
	void rotate(T const &angel, Vec3<T> const &axis)
	{
		m_stack.top() = m_stack.top() * Mat4<T>::rotate(angel,axis);
	}

	void push()
	{
		m_stack.push(m_stack.top());
	}
	void pop()
	{
		assert(m_stack.size()>1);
		m_stack.pop();
	}

	void triangle(Vec3<T> v[3])
	{

		for(int i=0;i<3;i++)
		{
			v[i] = transformPoint(v[i]);
			(*stage1)<<v[i]<<"\n";
			v[i] = m_view * v[i];
			(*stage2)<<v[i]<<"\n";
			v[i] = m_proj * v[i];
			(*stage3)<<v[i]<<"\n";
		}
		(*stage1)<<"\n";
		(*stage2)<<"\n";
		(*stage3)<<"\n";
		m_triangles.push_back(Triangle<T>(v,Color::random()));
	}

	Vec3<T> transformPoint(Vec3<T> p)
	{
		return m_stack.top()*p;
	}

	void setTriangles(ifstream & is)
	{
		m_triangles.clear();
		Vec3<T> v[3];
		while(is>>v[0])
		{
			for(int i=1;i<3;i++)
				is>>v[i];
			m_triangles.push_back(Triangle<T>(v,Color::random()));
		}
	}

	void draw_method1(
		int screen_width, 
		int screen_height,
		Vec3<T> mn,
		Vec3<T> mx,
		T dx,
		T dy,
		T top_y,
		T left_x,
		vector<vector<T> > &z_values,
		bitmap_image & image,
		Triangle<T> const& t
		)
	{
		// calculating normal of triangle
		Vec3<T> p = t[0];
		Vec3<T> n = (t[1]-t[0]).cross(t[2]-t[0]);

		// calculating coefficient of equation ax+by+cz+d=0;
		T a = n[0], b = n[1], c = n[2], d = -n.dot(p);

		// t.min(1) <=  top_y - i*dy <= t.max(1)
		int top_scanline = ceil((top_y - t.max(1))/dy);
		int bottom_scanline = floor((top_y - t.min(1))/dy);
		
		top_scanline = max(0,top_scanline);
		bottom_scanline = min(screen_height-1,bottom_scanline);

		// DBG(t.max(1));
		// DBG(t.min(1));
		// DBG(top_scanline);
		// DBG(bottom_scanline);
		
		T y_val = top_y - top_scanline*dy;
		for(
			int row = top_scanline;
			row<=bottom_scanline;
			row++,y_val -=dy
			)
		{
			
			assert(y_val<=t.max(1));
			assert(t.min(1)<=y_val);
			// if there is a 3d point inside triangle with y = y_val and minimize/maximize value of x
			// then that point should lie on border of triangles

			vector<T> x_values;
			Vec3<T> point = t[0];
			Vec3<T> vec = t[1]-t[0];
			{
				// vec.normalize();
				// point + a * vec  lies inside segment
				// now point.y + a * vec.y == y_val
				T a = (y_val - point[1])/vec[1];
				if(a>=T(0) and a<=T(1))
					x_values.push_back(point[0]+a*vec[0]);
			}

			vec = t[2]-t[0];
			{
				// vec.normalize();
				T a = (y_val - point[1])/vec[1];
				if(a>=T(0) and a<=T(1))
					x_values.push_back(point[0]+a*vec[0]);
			}
			
			point = t[1];
			vec = t[2]-t[1];
			{
				// vec.normalize();
				T a = (y_val - point[1])/vec[1];
				if(a>=T(0) and a<=T(1))
					x_values.push_back(point[0]+a*vec[0]);
			}
			assert(x_values.size());
			if(x_values.empty()) continue;

			// DBG(x_values.size());
			// for(auto i: x_values)
			// {
			// 	cerr<<i<<" ";
			// }
			// NL;
			
			sort(x_values.begin(),x_values.end());
			// DBG(dx);
			// DBG(dy);

			// x_min <= left_x + i * dx <= x_max
			int left_intersecting_col = ceil((*x_values.begin() - left_x)/dx);
			int right_intersecting_col = floor((*x_values.rbegin() - left_x)/dx); // dont use end() for last value
			left_intersecting_col = max(0,left_intersecting_col);
			right_intersecting_col = min(screen_width-1,right_intersecting_col);
			
			// DBG(left_intersecting_col);
			// DBG(right_intersecting_col);

			T x_val = left_x + left_intersecting_col*dx;
			T z_value = - (a*x_val + d+b*y_val)/c;
			T z_incr = -a*dx/c;
			for(
				int col=left_intersecting_col;
				col<=right_intersecting_col ; 
				col++,z_value += z_incr
				,x_val += dx
				)
			{
				//checking if x_value are in range
				assert(x_val<=*x_values.rbegin());
				assert(*x_values.begin()<=x_val);

				// Calculate z value from triangle t
				if(z_values[col][row] <= z_value) continue;
				if(z_value<mn[2] ) continue;
				image.set_pixel(col,row,t.c[0],t.c[1],t.c[2]);
				z_values[col][row] = z_value;
			}
		}
	}
	


	void draw_method2(
		int screen_width, 
		int screen_height,
		Vec3<T> mn,
		Vec3<T> mx,
		T dx,
		T dy,
		T top_y,
		T left_x,
		vector<vector<T> > &z_values,
		bitmap_image & image,
		Triangle<T> const& t
		)
	{
		// y = top_y - i*dy , 0<=i < screen_height
		// x = left_x + i*dx , 0<=i < screen_width
		int top_scanline = max(0,(int)(ceil((top_y-t.max(1))/dy)) );
		int bottom_scanline = min(screen_height-1, (int)(floor((top_y-t.min(1))/dy)));
		T ys = top_y - top_scanline*dy;
		for(
			int row = top_scanline;
			row<=bottom_scanline;
			row++, ys-=dy
			)
		{
			//      t[0]
			//      / \
			//     /   \
			//    /     \
			// t[1]      \
			//    .       \
			//      .      \
			//         .    \
			//           . t[2]
			// 
			vector<pair<T,int> > x_values;
			for(int i=0;i<3;i++)
			{
				// checking intersection point on segment [t[i],t[i+1])
				
				int i_1 = (i+1)%3;
				if(t[i][1] == t[i_1][1]) continue;
				if(ys < min(t[i][1],t[i_1][1])) continue;
				if(max(t[i][1],t[i_1][1]) < ys) continue;

				T x = t[i][0] + (ys - t[i][1])*(t[i_1][0]-t[i][0])/(t[i_1][1]-t[i][1]);
				if(x < min(t[i][0],t[i_1][0])) continue;
				if(max(t[i][0],t[i_1][0]) < x) continue;
				x_values.push_back(make_pair(x,i));
			}
			assert(!x_values.empty());

			assert(x_values.size()!=1);
			assert(x_values.size()!=3); // can happen 
			// line1           line2        intersection  other two
			// line[0 -> 1] & line[1 -> 2] -> p1 = 1 , (p2,p3) = (0,2)
			// line[0 -> 1] & line[2 -> 0] -> p1 = 0 , (p2,p3) = (1,2)
			// line[1 -> 2] & line[2 -> 0] -> p1 = 2 , (p2,p3) = (0,1)
			int intersec = -1;
			// x_values are sorted by line index[second]
			if(x_values[0].second == 0)
				intersec = x_values[1].second == 1? 1:0;
			else
				intersec = 2;
			int other_1 = (intersec+1)%3;
			int other_2 = (intersec+2)%3;

			sort(x_values.begin(),x_values.end());

			// left_x + i * dy = x
			int left_intersecting_column = max(
					0,
					(int) ceil( (x_values.begin()->first - left_x)/dx )
				);
			int right_intersecting_column = min(
					screen_width-1,
					(int) floor( (x_values.rbegin()->first - left_x)/dx )
				);
			T x_a = t[intersec][0] + (ys-t[intersec][1])*
								(t[other_1][0]-t[intersec][0])/
								(t[other_1][1]-t[intersec][1]);
			T x_b = t[intersec][0] + (ys-t[intersec][1])*
								(t[other_2][0]-t[intersec][0])/
								(t[other_2][1]-t[intersec][1]);
			T z_a = t[intersec][2] + (ys-t[intersec][1])*
								(t[other_1][2]-t[intersec][2])/
								(t[other_1][1]-t[intersec][1]);
			T z_b = t[intersec][2] + (ys-t[intersec][1])*
								(t[other_2][2]-t[intersec][2])/
								(t[other_2][1]-t[intersec][1]);
			if(x_a > x_b)
			{
				swap(x_a,x_b);
				swap(z_a,z_b);
				swap(other_1,other_2);
			}
			assert(x_a <= x_b);
			assert(abs(x_a - x_values.begin()->first)<=EPS);
			assert(abs(x_b - x_values.rbegin()->first)<=EPS);
			T x_p = left_x + left_intersecting_column*dx;
			T z_p = z_a + (x_p - x_a)*(z_b-z_a)/(x_b-x_a);
			T dz = dx*(z_b-z_a)/(x_b-x_a);
			for(
				int column = left_intersecting_column;
				column<=right_intersecting_column;
				column++,
				x_p+=dx,
				z_p+=dz
				)
			{
				assert(x_p >= x_a);
				assert(x_p <= x_b);
				if(z_values[column][row] <= z_p) continue;
				if(z_p < mn[2]) continue;
				z_values[column][row] = z_p;
				image.set_pixel(column,row,t.c[0],t.c[1],t.c[2]);
			}
		}

	}
	
	void draw(int screen_width, 
				int screen_height,
				Vec3<T> mn,
				Vec3<T> mx,
				string output_file_name
				)
	{
		
		T dx = (mx[0]-mn[0])/T(screen_width);
		T dy = (mx[1]-mn[1])/T(screen_height);
		
		T top_y = mx[1] - dy/2;
		T left_x = mn[0] + dx/2;

		// init z-buffer
		vector< vector<T> > z_values(screen_width,vector<T>(screen_height,mx[2]));
		bitmap_image image(screen_width,screen_height);

		for(int i=0;i<screen_width;i++)
		{
			for(int j=0;j<screen_height;j++)
			{
				// init color as black
				image.set_pixel(i,j,0,0,0);
			}
		}

		// DBG(top_y);
		// DBG(left_x);
		// DBG(dy);

		for(Triangle<T> const&t:m_triangles)
		{
			draw_method2(
				screen_width,
				screen_height,
				mn,
				mx,
				dx,
				dy,
				top_y,
				left_x,
				z_values,
				image,
				t
				);
		}
		
		for(int i=0;i<screen_height;i++)
		{
			for(int j=0;j<screen_width;j++)
			{
				if(z_values[j][i] < mx[2])
					(*z_out)<<z_values[j][i]<<"\t";
				// else 
				// 	(*z_out)<<string(z_out->precision(),' ')<<" ";
			}
			(*z_out)<<"\n";
		}

		image.save_image(output_file_name);
	}
	
};

class Color
{
	private:
		int c[3];
	public:
		Color(int r=0,int g=0,int b=0):c{r,g,b}{
			
		}
		static Color red(){
			return Color(255,0,0);
		}
		static Color green(){
			return Color(0,255,0);
		}
		static Color blue(){
			return Color(0,0,255);
		}
		static Color black(){
			return Color(0,0,0);
		}
		static Color white(){
			return Color(255,255,255);
		}
		static Color random(){
			return Color(rand()%255,rand()%255,rand()%255);
		}

		// << operator overloading
		friend ostream &operator <<(ostream &os , Color &c)
		{
			os<<"("<<c.c[0]<<" "<<c.c[1]<<" "<<c.c[2]<<")";
			return os;
		}

		// [] operator overloading
		int & operator[](int i) const
		{
			return const_cast< int &>(c[i]);
		}

};
template<class T>
class Triangle
{
	Vec3<T> v[3];
public:
	Color c;
	Triangle(Vec3<T> v[3],Color c):c(c){
		for(int i=0;i<3;i++){
			this->v[i] = v[i];
		}
	}
	T min(int i) const
	{
		T ret = v[0][i];
		if(v[1][i]<ret)
			ret = v[1][i];
		if(v[2][i]<ret)
			ret = v[2][i];
		return ret;
	}

	T max(int i) const
	{
		T ret = v[0][i];
		if(v[1][i]>ret)
			ret = v[1][i];
		if(v[2][i]>ret)
			ret = v[2][i];
		return ret;
	}

	Vec3<T> & operator[](int i) const
	{
		return const_cast< Vec3<T> &>(v[i]);
	}
	
	friend ostream &operator <<(ostream &os , Triangle<T> &t)
	{
		os<<"Triangle: "<<"\n";
		for(int i=0;i<3;i++)
			os<<t.v[i]<<"\n";
		os<<"Color: "<<t.c<<"\n";
		return os;
	}

};


template<typename T>
class Vec4
{
	T v[4];
public:

	// Vec4(T x=T(0), T y=T(0), T z=T(0), T w=T(1))
	Vec4(T x, T y, T z, T w)
	{
		v[0] = x;
		v[1] = y;
		v[2] = z;
		v[3] = w;
	}
	T& operator[](int const &x) const
	{
		return const_cast<T&>(v[x]);
	}
	friend ostream &operator<<(ostream &os, Vec4<T> const &v4)
	{
		os << v4.v[0] << " " << v4.v[1] << " " << v4.v[2] << " " << v4.v[3];
		return os;
	}
	
};


template<typename T>
class Vec3: public Vec4<T> 
{
	public:
	Vec3(T x=T(0), T y=T(0), T z=T(0)):Vec4<T>(x, y, z, T(1))
	{
		
	}
	T length() const
	{
		T len = T(0);
		for(int i=0;i<3;i++)
		{
			len += this->operator[](i) * this->operator[](i);
		}
		len = sqrt(len);
		return len;
	}
	void normalize()
	{
		T len = length();
		for(int i=0;i<3;i++)
		{
			this->operator[](i) /= len;
		}
	}
	T dot(Vec3<T> const & v) const
	{
		T ret = T(0);
		for(int i=0;i<3;i++)
		{
			ret += this->operator[](i) * v[i];
		}
		return ret;
	}
	Vec3<T> cross(Vec3<T> const & v) const
	{
		
		Vec3<T> ret; // (0,0,0,1)
		#ifdef X
			#error "X is defined"
		#endif

		#define X 0
		#define Y 1
		#define Z 2

		ret[X] = this->operator[](Y)*v[Z] - this->operator[](Z)*v[Y];
		ret[Y] = this->operator[](Z)*v[X] - this->operator[](X)*v[Z];
		ret[Z] = this->operator[](X)*v[Y] - this->operator[](Y)*v[X];

		#undef X
		#undef Y
		#undef Z
		return ret;
	}
	Vec3<T> rotate(T const& angel,Vec3<T> const & axis)
	{
		Vec3<T> ret; // (0,0,0,1)
		T cos_theta = cos(angel* PI / 180.0);
		ret += (*this * cos_theta);
		ret += axis * (axis.dot(*this)*(1-cos_theta));
		ret += axis.cross(*this) * sin(angel* PI / 180.0);
		
		return ret;	
	}

	Vec3<T> operator - (Vec3<T> const & v) const
	{
		Vec3<T> ret; // (0,0,0,1)
		for(int i=0;i<3;i++)
		{
			ret[i] = this->operator[](i) - v[i];
		}
		return ret;
	}
	Vec3<T> operator-()const
	{
		Vec3<T> ret; // (0,0,0,1)
		for(int i=0;i<3;i++)
		{
			ret[i] = -this->operator[](i);
		}
		return ret;
	}

	Vec3<T> operator *(T const & m) const
	{
		Vec3<T> ret; // (0,0,0,1)
		for(int i=0;i<3;i++)
		{
			ret[i] = this->operator[](i) * m;
		}
		return ret;
	}
	Vec3<T> &operator+=(Vec3<T> const & v)
	{
		for(int i=0;i<3;i++)
		{
			this->operator[](i) += v[i];
		}
		return *this;
	}
	
	T& operator[](int const & x) const
	{
		return Vec4<T>::operator[](x);
	}

	friend istream & operator >>(istream & is,Vec3<T> & v3)
	{
		is>>v3[0]>>v3[1]>>v3[2];
		v3[3]=T(1);
		return is;
	}

	friend ostream & operator <<(ostream & os,Vec3<T> const & v3)
	{
		return os<<v3[0]<<" "<<v3[1]<<" "<<v3[2];
	}
};

template<typename T>
class Vec4
{
	T v[4];
public:

	// Vec4(T x=T(0), T y=T(0), T z=T(0), T w=T(1))
	Vec4(T x, T y, T z, T w)
	{
		v[0] = x;
		v[1] = y;
		v[2] = z;
		v[3] = w;
	}
	T& operator[](int const &x) const
	{
		return const_cast<T&>(v[x]);
	}
	friend ostream &operator<<(ostream &os, Vec4<T> const &v4)
	{
		os << v4.v[0] << " " << v4.v[1] << " " << v4.v[2] << " " << v4.v[3];
		return os;
	}
	
};
template<typename T>
class Mat4
{
	T m[4][4];
	public:
	/*
		Zero matrix
	*/
	Mat4()
	{
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				m[i][j] = T(0);
			}
		}
	}
	/*
		Identity matrix
	*/
	static Mat4 identity()
	{
		Mat4<T> ret; // zero matrix
		for(int i=0;i<4;i++)
		{
			ret[i][i]=1;
		}
		return ret;
	}
	
	static Mat4 translate(Vec3<T> const& v)
	{
		Mat4<T> ret = identity();
		ret[0][3] = v[0];
		ret[1][3] = v[1];
		ret[2][3] = v[2];
		return ret;
	}

	static Mat4 scale(Vec3<T> &v)
	{
		Mat4<T> ret ; // zero
		ret[0][0] = v[0];
		ret[1][1] = v[1];
		ret[2][2] = v[2];
		ret[3][3] = T(1);
		return ret;
	}

	static Mat4 rotate(T const & angel, Vec3<T> const & axis)
	{
		Mat4<T> ret; // zero
		Vec3<T> a = axis;
		a.normalize();

		ret[3][3] = T(1);
		for(int i=0;i<3;i++)
		{
			Vec3<T> unit; // (0,0,0,1)
			unit[i]=T(1);
			Vec3<T> rotated = unit.rotate(angel,a);
			for(int j=0;j<3;j++)
				ret[j][i] = rotated[j];
		}
		return ret;
	}


	
	T * operator [](int const & x) const
	{
		return const_cast<T*>(m[x]);
	}

	Mat4 operator*(Mat4 const& rhs) const
	{
		Mat4 ret; // zero matrix
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				for(int k=0;k<4;k++)
				{
					ret[i][j] += m[i][k]*rhs[k][j];
				}
			}
		}
		return ret;
	}
	
	Vec3<T> operator *(Vec3<T> const& vec3) const
	{
		Vec3<T> ret; // (0,0,0,1)
		ret[3]=0;
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				ret[i] += m[i][j]*(vec3[j]);
			}
		}
		for(int i=0;i<3;i++)
			ret[i]/=ret[3];
		ret[3]=T(1);
		return ret;
	}



	void to_3x3(ostream &os)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				os<<m[i][j]<<" ";
			}
			os<<"\n";
		}
	}

	friend ostream & operator <<(ostream & os,Mat4 const & m4)
	{
		os<<"\n";
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++)
				os<<m4.m[i][j]<<" ";
			os<<"\n";
		}
		return os;
	}


};

#endif // GL_H
