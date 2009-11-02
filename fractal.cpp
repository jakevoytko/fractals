#include <cmath>
#include <cstdlib>
#include <deque>
#include <GL/glut.h>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

#define window_width  640
#define window_height 480

const double pi = 3.141592653589793;

enum Fractal {
    FractalBoxOutline,
    FractalDragonCurve,
    FractalExteriorSnowflake,
    FractalHilbertICurve,
    FractalHilbertIICurve,
    FractalKochAntisnowflake,
    FractalKochCurve,
    FractalKochIsland,
    FractalKochSnowflake,
    FractalPeanoCurve,
    FractalSierpinskiArrowhead,
    FractalSierpinskiCurve,
    FractalSierpinskiTriangle,
    FractalNumFractals}; // Not an actual fractal, used for % op

// My kingdom for inline dictionary syntax
string FractalName(Fractal f)
{
    switch(f){
    case FractalBoxOutline:
	return "Box outline";
    case FractalDragonCurve:
	return "Dragon curve";
    case FractalExteriorSnowflake:
	return "Exterior Snowflake";
    case FractalHilbertICurve:
	return "Hilbert I Curve";
    case FractalHilbertIICurve:
	return "Hilbert II Curve";
    case FractalKochAntisnowflake:
	return "Koch Antisnowflake";
    case FractalKochCurve:
	return "Koch curve";
    case FractalKochIsland:
	return "Koch island";
    case FractalPeanoCurve:
	return "Peano curve";
    case FractalSierpinskiArrowhead:
	return "Sierpinski Arrowhead";
    case FractalSierpinskiCurve:
	return "Sierpinski Curve";
    case FractalSierpinskiTriangle:
	return "Sierpinski Triangle";
    case FractalKochSnowflake:
	return "Koch Snowflake";
    default:
	return "";
    }
}

Fractal fractal = FractalKochSnowflake;

int level = 3;

template <class T>
class vertex2
{
public:
    vertex2(T x, T y)
	:m_x(x), m_y(y) {
    }
    
    T x() const {return m_x;}
    T y() const {return m_y;}
private:
    T m_x;
    T m_y;
};


template <class T>
class vect2 : public vertex2<T>
{
public:
    vect2(T x, T y)
	:vertex2<T>(x, y){
    }
    
    T angle() const { return atan(this->y() / this->x()); }
    T magnitude() const {return sqrt(this->x() * this->x()
				     + this->y() * this->y()); }
};


template <class T>
class ray2
{
public:
    ray2(const vect2<T> &vect, const vertex2<T> &anchor)
	:m_vect(vect), m_anchor(anchor) {
    }

    vect2<T> vect() const { return m_vect; }
    vertex2<T> anchor() const { return m_anchor; }
    
private:
    vect2<T> m_vect;
    vertex2<T> m_anchor;
};


template <class T>
vertex2<T> operator+(const vertex2<T> &v, const vect2<T> &vect)
{
    return vertex2<T>(v.x() + vect.x(), v.y() + vect.y());
}

template <class T>
ray2<T> operator+(const ray2<T> &ray, const vect2<T> &vect)
{
    return ray2<T>(ray.vect(), ray.anchor() + vect);
}

template <class T>
vect2<T> operator*(const vect2<T> &vect, T constant)
{
    return vect2<T>(vect.x() * constant, vect.y() * constant);
}

template <class T>
ray2<T> operator*(const ray2<T> &ray, T constant)
{
    return ray2<T>(ray.vect() * constant, ray.anchor());
}

// Returns a new ray2 where the anchor equals the sum of the old anchor
// and the old vector. The new vector equals the old vector.
template <class T>
ray2<T> Advance(const ray2<T> &ray)
{
    return ray2<T>(ray.vect(), ray.anchor() + ray.vect());
}

typedef vertex2<float> vertex2f;
typedef ray2<float> ray2f;

void Text(const vertex2<float> &pt, const string& text)
{
    bool blend = glIsEnabled(GL_BLEND);
    glEnable(GL_BLEND);

    glRasterPos2f(pt.x(), pt.y());

    for(int i=0; i<text.size(); ++i)
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);

    if(!blend)
	glDisable(GL_BLEND);
}

void Line(const vertex2<float> &pt1, const vertex2<float> &pt2)
{
    glBegin(GL_LINES);
    glVertex2f((float)pt1.x(), (float)pt1.y());
    glVertex2f((float)pt2.x(), (float)pt2.y());
    glEnd();
}

void Line(const ray2<float> &ray)
{
    vertex2<float> anchor = ray.anchor();
    vertex2<float> dest = anchor + ray.vect();
    
    glBegin(GL_LINES);
    glVertex2f(anchor.x(), anchor.y());
    glVertex2f(dest.x(), dest.y());
    glEnd();
}

template <class T>
vect2<T> RotateCw(const vect2<T> &vect, T radians)
{
    T x = cos(radians)*vect.x() + sin(radians)*vect.y();
    T y = -sin(radians)*vect.x() + cos(radians)*vect.y();
    return vect2<T>(x, y);
}

template <class T>
vect2<T> RotateCcw(const vect2<float> &vect, T radians)
{
    return RotateCw(vect, -radians);
}

template <class T>
ray2<T> RotateCw(const ray2<T> &ray, T radians)
{
    return ray2<T>(RotateCw(ray.vect(), radians), ray.anchor());
}

template <class T>
ray2<T> RotateCcw(const ray2<float> &ray, T radians)
{
    return ray2<T>(RotateCcw(ray.vect(), radians), ray.anchor());
}

//
// Stores information about a Lindenmayer system.
// Call LSystem::Level() to get the expanded system up to a certain
// fractal level
//
// In any given representation string, the rules contain terminals '+'
// and '-'. '+' means turn CW, '-' means turn CCW.  Nonterminal 'F'
// means 'move forward'. All other symbols should be ignored for
// walks.
class LSystem
{
public:
    LSystem() // ARGH required by std::map. Honestly, C++.
    {
	m_levels.push_back("F++F++F");
	m_rules['F'] = "F-F++F-F";
	m_minusTurn = m_plusTurn = (float)(pi/3);
    }

    // level0: The string representing the base case
    // rules: Translation rules.
    LSystem(const string& level0, const map<char, string> rules,
	    float minusTurn, float plusTurn, float K = (float)(1./3.)):
	m_rules(rules),m_minusTurn(minusTurn), m_plusTurn(plusTurn),
	m_K(K)
    {
	m_levels.push_back(level0);
    }

    // Returns the representation after 'level' string applications.
    std::string Level(int level){
	if(level > m_levels.size() - 1){
	    for(int i=m_levels.size(); i<=level; ++i){
		m_levels.push_back(GenerateNextLevel(m_levels[i-1]));
	    }
	}
	return m_levels[level];
    }

    // The angle to turn CW
    float MinusAngle() { return m_minusTurn; }

    // The angle to turn CCW
    float PlusAngle() { return m_plusTurn; }

    // The constant by which 'F' shrinks per number of iterations.
    // Sierpinski trangle = .5, snowflake = 1/3, etc
    float K() { return m_K; }
private:

    std::string GenerateNextLevel(const string &tape)
    {
	ostringstream fmt;
	for(int i=0; i<tape.size(); ++i){
	    if(m_rules.find(tape[i]) != m_rules.end())
		fmt<<m_rules[tape[i]];
	    else
		fmt<<tape[i];
	}
	return fmt.str();
    }

    deque<string> m_levels;
    map<char, string> m_rules;
    float m_minusTurn, m_plusTurn;
    float m_K;
};

map<Fractal, LSystem> fractalMap;

void WalkLSystem(const string &lSystem, const ray2f F,
		 float minusAngle, float plusAngle)
{
    ray2f current = F;
    for(int i=0; i<lSystem.size(); ++i){
	switch( lSystem[i]){
	case '-':
	    current = RotateCcw(current, minusAngle);
	    break;
	case '+':
	    current = RotateCw(current, plusAngle);
	    break;
	case 'F':
	    Line(current);
	    current = Advance(current);
	    break;
	}
    }
}

LSystem KochSnowflake()
{
    map<char, string> rules;
    rules['F'] = "F-F++F-F";
    return LSystem("F++F++F", rules, (float)(pi/3.), (float)(pi/3.));
}

LSystem KochAntisnowflake()
{
    map<char, string> rules;
    rules['F'] = "F+F-F+F";
    return LSystem("F++F++F++F", rules, (float)(2*pi/3), (float)(pi/3));
}


LSystem BoxOutline()
{
    map<char, string> rules;
    rules['F'] = "F+F-F-F+F";
    return LSystem("F+F+F+F", rules, (float)(pi/2), (float)(pi/2));
}

LSystem SierpinskiTriangle()
{
    map<char, string> rules;
    rules['F'] = "FF";
    rules['X'] = "++FXF--FXF--FXF++";
    return LSystem("FXF++FF++FF", rules, (float)(pi/3.), (float)(pi/3.), .5);
}

LSystem KochCurve()
{
    map<char, string> rules;
    rules['F'] = "F+F-F-F+F";
    return LSystem("F", rules, (float)(pi/2), (float)(pi/2));
}

LSystem KochIsland()
{
    map<char, string> rules;
    rules['F'] = "F+F-F-FF+F+F-F";
    return LSystem("F+F+F+F", rules, (float)(pi/2), (float)(pi/2), .25);
}

LSystem SierpinskiArrowhead()
{
    map<char, string> rules;
    rules['X'] = "YF+XF+Y";
    rules['Y'] = "XF-YF-X";
    return LSystem("YF", rules, (float)(pi/3), (float)(pi/3), .5);
}

LSystem DragonCurve()
{
    map<char, string> rules;
    rules['X'] = "X+YF+";
    rules['Y'] = "-FX-Y";
    return LSystem("YF", rules, (float)(pi/2), (float)(pi/2), .7);
}

LSystem HilbertICurve()
{
    map<char, string> rules;
    rules['L'] = "+RF-LFL-FR+";
    rules['R'] = "-LF+RFR+FL-";
    return LSystem("L", rules, (float)(pi/2), (float)(pi/2), .5);
}

LSystem HilbertIICurve()
{
    map<char, string> rules;
    rules['X'] = "XFYFX+F+YFXFY-F-XFYFX";
    rules['Y'] = "YFXFY-F-XFYFX+F+YFXFY";
    return LSystem("X", rules, (float)(pi/2), (float)(pi/2));
}

LSystem ExteriorSnowflake()
{
    map<char, string> rules;
    rules['F'] = "F+F-F+F";
    return LSystem("F+F+F+F+F+F", rules, (float)(2*pi/3), (float)(pi/3));
}

LSystem PeanoCurve()
{
    map<char, string> rules;
    rules['F'] = "F-F+F+FF+F+F+FF";
    return LSystem("F", rules, (float)(pi/2), (float)(pi/2));
}

LSystem SierpinskiCurve()
{
    map<char, string> rules;
    rules['F'] = "F+F-F+F-F";
    return LSystem("F+F+F+F", rules, (float)(pi/2), (float)(pi/2), .25);
}

// Main loop
void main_loop_function()
{
    // Initialize screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0,0, -5);

    ray2f walk(vect2<float>(1.4, 0), vertex2<float>(-.5, .5));

    string walkRule = fractalMap[fractal].Level(level);
    float minusAngle = fractalMap[fractal].MinusAngle();
    float plusAngle = fractalMap[fractal].PlusAngle();
    float K = fractalMap[fractal].K();
    glColor4f(0, 1, 0, 1);
    Text(vertex2<float>(0, -2), FractalName(fractal));
    Text(vertex2<float>(0, 1.9), string("Level ") + (char)(((int)level)+'0'));
    glColor4f(0, 0, 1, 1);
    WalkLSystem(walkRule, walk * (float)(pow(K, level)),
		minusAngle, plusAngle);

    glutSwapBuffers();
}

// Initialze OpenGL perspective matrix
void GL_Setup(int width, int height)
{
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glEnable( GL_DEPTH_TEST );
    gluPerspective( 45, (float)width/height, .1, 100 );
    glMatrixMode( GL_MODELVIEW );
}

void processNormalKeys(unsigned char key, int x, int y) {
    int numFractals = (int)FractalNumFractals;
    int currentFractal = (int)fractal;
    
    if (key == 'q' || key == 'Q') 
	exit(0);
    else if(key >= '0' && key <= '9')
	level = key - '0';
    // Playing some hackey sack to make this easy
    else if(key == '+' || key == '=')
	fractal = ((Fractal)((currentFractal+1) % numFractals));
    else if(key == '-' || key == '_'){ 
	// Why doesn't % act like it does in number theory?!
	--currentFractal;
	if(currentFractal < 0)
	    currentFractal = numFractals - 1;
	fractal = ((Fractal)(currentFractal));
    }
}

void CreateLSystems()
{
    fractalMap[FractalBoxOutline] = BoxOutline();
    fractalMap[FractalDragonCurve] = DragonCurve();
    fractalMap[FractalExteriorSnowflake] = ExteriorSnowflake();
    fractalMap[FractalHilbertICurve] = HilbertICurve();
    fractalMap[FractalHilbertIICurve] = HilbertIICurve();
    fractalMap[FractalKochCurve] = KochCurve();
    fractalMap[FractalKochAntisnowflake] = KochAntisnowflake();
    fractalMap[FractalKochIsland] = KochIsland();
    fractalMap[FractalKochSnowflake] = KochSnowflake();
    fractalMap[FractalPeanoCurve] = PeanoCurve();
    fractalMap[FractalSierpinskiArrowhead] = SierpinskiArrowhead();
    fractalMap[FractalSierpinskiCurve] = SierpinskiCurve();
    fractalMap[FractalSierpinskiTriangle] = SierpinskiTriangle();
}


// Initialize GLUT and start main loop
int main(int argc, char** argv) {
    CreateLSystems();
    glutInit(&argc, argv);
    glutInitWindowSize(window_width, window_height);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutCreateWindow("Fractal demo");
    glutKeyboardFunc(processNormalKeys);
    glutIdleFunc(main_loop_function);
    GL_Setup(window_width, window_height);
    glutMainLoop();
}
