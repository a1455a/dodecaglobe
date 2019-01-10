//
// Read world map in Miller cylindrical progection from SVG,
// write one SVG projection for paper dodecaglobe
//
#include <expat.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

// Output size and formatting
const double paper_width_in_inches = 8.5;
const double paper_height_in_inches = 11;
const double paper_margin_in_inches = 1;
const double pentagon_radius_in_inches = 2.1;
const int dpi = 100;

// Determined by trial and error to approximately match Google Maps
const double x_zero = 315.91904095;
const double y_zero = 500.59241;
const double target_width = 3389.83;
const double target_height = 2536.55;

const double parallel49n = 715; // border between US and Canada
const double meridian141w = 274;

const double meridian141e = 2925.5; // border between Indonesia and Papua New Guinea

const double greenwich = (meridian141w + meridian141e) / 2;
const double equator = 1255;

const double parallel22s = 1484; // border between Namibia and Botswana
const double meridian20e = greenwich + (meridian141e - greenwich) * 20 / 141;

const double parallel22n = 1025.5 ; // southwest corner of Egypt
const double meridian25e = greenwich + (meridian141e - greenwich) * 25 / 141;

// Dimensions of dodecahedron face inscribed in unit sphere
const double h = sqrt((5 + sqrt(5)) / 30); // Distance between face center and edge ("feet")
const double r = sqrt(2 * (5 - sqrt(5)) / 15); // Distance between face center and vertex ("head")
const double s = sqrt(2 * (3 - sqrt(5)) / 3); // Side
const double w = 2 * s * sin(54 * M_PI / 180); // Width

// Often used calculated values
const double printable_width_in_px = (paper_width_in_inches - 2 * paper_margin_in_inches) * dpi;
const double printable_height_in_px = (paper_height_in_inches - 2 * paper_margin_in_inches) * dpi;
const double pentagon_radius_in_pixels = pentagon_radius_in_inches * dpi;
const double output_scaling_factor = pentagon_radius_in_pixels / r;
const std::complex<double> pentagon_center( printable_width_in_px / 2 + h * output_scaling_factor * sin(18 * M_PI / 180)
	                                      , printable_height_in_px / 2 - h * output_scaling_factor * cos(18 * M_PI / 180));
const std::complex<double> mirror_center( printable_width_in_px / 2 - h * output_scaling_factor * sin(18 * M_PI / 180)
	                                    , printable_height_in_px / 2 + h * output_scaling_factor * cos(18 * M_PI / 180));
 
// In a small stand-alone program we can get away with using globals instead of passing around struct UserData*
double lambda0; // Rotate the Earth to align zero medidian with one of the edges
double a = 1, b, c, d = 1, e, f; // Some shapes have "transform=matrix(a,b,c,d,e,f)"
char country_cd_1, country_cd_2; // Tho-character country code
bool south, polar, land, water;
double face_rotation_angle = 0; // Orientation of face on paper

// Convert pixel coordinates of the flat map to pixel coordinates of a printed face.
// Return x coordinate in real part and y coordinate in imaginary part of a complex number.
std::complex<double> project(double x, double y) {

  // Convert to unitless x and y, then undo Miller cylindrical projection to get 
  // lambda (longitude) and phi (latitude), in radians
  double lambda = M_PI * (x - greenwich) * (141 + 141) / ((meridian141e - meridian141w) * 180) - lambda0;
  const double yfactor = (equator - parallel49n) / (5 * log(tan(M_PI / 4 + 2 * 49 * M_PI / (180 * 5)))/ 4);
  double phi_sign = ((y < equator) ^ south ? 1 : -1);

  double phi = phi_sign * (5 * atan(exp(4 * fabs(y - equator) / (5 * yfactor))) / 2 - 5 * M_PI / 8);

  // Dihedral angle of regular dodecahedron (angle between adjacent faces)
  const double psi = acos(-1 / sqrt(5));

  // Distance from the center of the regular dodecahedron inscribed in unit sphere to a face
  const double p = sqrt((5 + sqrt(5)) / 30) * tan(psi / 2);

  // Distance of projected point from the face center
  double projected_x = nan(""), projected_y = nan("");
  if (polar) {
  	if (phi > M_PI / 4) {
      projected_x = -p * sin(lambda) / tan(phi);
      projected_y = (south ? 1 : -1) * p * cos(lambda) / tan(phi);
    }
  } else {
    if ( lambda > -M_PI / 2 && lambda < M_PI / 2
  	  && phi > -M_PI / 4 && phi < 7 * M_PI / 8
  	   ) {
      double eta = atan(tan(phi) / cos(lambda));
      projected_x = (south ? -1 : 1) * p  * tan(lambda) * cos(eta) / cos(eta - psi + M_PI / 2);
      projected_y = -p * tan(eta - psi + M_PI / 2);
    }
  }

  // Rotate and scale
  const std::complex<double> rotation(cos(face_rotation_angle), -sin(face_rotation_angle));
  const std::complex<double> projected(projected_x, projected_y);
  return (pentagon_center + projected * rotation * output_scaling_factor);
}

static void XMLCALL startElement(void* userData, const char* name, const char** attr)
{
  const char** p = attr;
  if (!strcmp(name, "g")) {
    while (*p) {
      if (!strcmp(*p, "class")) {
      	if (!strcmp(p[1], "ocean")) {
          land = false;
          water = true;
      	} else if (!strcmp(p[1], "lake")) {
          land = false;
          water = true;
      	} else if (p[1][0] == 'c') {
      	  land = false;
      	  water = false;
      	} else {
      	  if (strlen(p[1]) == 7) { // "land us"
      	  	country_cd_1 = p[1][5];
      	  	country_cd_2 = p[1][6];
      	  }
      	  land = true;
      	  water = false;
      	}
      } else if (!strcmp(*p, "transform")) {
        sscanf(p[1], "matrix(%lf,%lf,%lf,%lf,%lf,%lf)", &a, &b, &c, &d, &e, &f);
      }
      p += 2;
    }
  } else if (!strcmp(name, "path")) {
    while (*p) {
      if (!strcmp(*p, "class")) {
      	if (!strcmp(p[1], "ocean")) {
          land = false;
          water = true;
      	} else if (!strcmp(p[1], "lake")) {
          land = false;
          water = true;
      	} else if (p[1][0] == 'c') {
      	  land = false;
      	  water = false;
      	} else {
      	  if (strlen(p[1]) == 7) { // "land us"
      	  	country_cd_1 = p[1][5];
      	  	country_cd_2 = p[1][6];
      	  }
      	  land = true;
      	  water = false;
      	}
      } else if ((land || water) && !strcmp(*p, "d")) {

        const char* q = p[1]; // M x,y C x1,y1 x2,y2 x3,y3
        double sx, sy;
        if (sscanf(p[1], "M %lf,%lf", &sx, &sy) == 2) {

          // Apply transformation matrix
          double tsx = a * sx + c * sy + e + x_zero;
          double tsy = b * sx + d * sy + f + y_zero;

          // Project on dodecadron face
          std::complex<double> ts(project(tsx, tsy));

          // Special cases
          if (!std::isfinite(ts.real())) {
            // Russia starts on the other side of the world, use a closer starting point
            // so that Kamchatka is shown
          	if (country_cd_1 == 'r' 
          	  && country_cd_2 == 'u' 
          	  && lambda0 ==  4 * M_PI / 5 + M_PI / 5
          	  ) {
              ts = std::complex<double>(meridian141e, parallel49n);
            }
            // Same with Canada
          	if (country_cd_1 == 'c' 
          	  && country_cd_2 == 'a' 
          	  && lambda0 ==  -2 * M_PI / 5 + M_PI / 5
          	  ) {
              ts = std::complex<double>(meridian141w, parallel49n);
            }
          }

          if (std::isfinite(ts.real())) { // Ignore shapes from the other side of the globe
            double ax, ay, bx, by, cx, cy;
            bool may_be_visible = false;

            q = strchr(q, 'C');
            if (q) {

              std::vector<std::complex<double> > av, bv, cv;

              while (sscanf(q, "C %lf,%lf %lf,%lf %lf,%lf", &ax, &ay, &bx, &by, &cx, &cy) == 6) {

                // Apply transformation matrix
                double wax = a * ax + c * ay + e + x_zero;
                double way = b * ax + d * ay + f + y_zero;
                double wbx = a * bx + c * by + e + x_zero;
                double wby = b * bx + d * by + f + y_zero;
                double wcx = a * cx + c * cy + e + x_zero;
                double wcy = b * cx + d * cy + f + y_zero;

                // Project on dodecadron face
                std::complex<double> pa(project(wax, way)); 
                std::complex<double> pb(project(wbx, wby));
                std::complex<double> pc(project(wcx, wcy));

                if ( std::isfinite(pa.real()) 
                  && std::isfinite(pb.real())
                  && std::isfinite(pc.real())
                   ) {
                  av.push_back(pa);
                  bv.push_back(pb);
                  cv.push_back(pc);

                  if ( std::abs(pa - pentagon_center) < pentagon_radius_in_pixels
                  	|| std::abs(pc - pentagon_center) < pentagon_radius_in_pixels
                  	 ) {
                  	may_be_visible = true;
                  }
                }
                q = strchr(q + 1, 'C');

                if (!q) break;
              }
              if (may_be_visible) {
                std::cout << "<path ";
                if (country_cd_1 && country_cd_2) {
                  std::cout << "class=\"" << country_cd_1 << country_cd_2 << "\" ";
                }
                std::cout << "d=\"M " << ts.real() << "," << ts.imag();
                for (int i = 0; i < av.size(); i++) {
                   std::cout << " C " << av[i].real() << "," << av[i].imag() << " " 
                                      << bv[i].real() << "," << bv[i].imag() << " " 
                                      << cv[i].real() << "," << cv[i].imag();
                }
                std::cout << "\" fill=\"#" << (water ? "ffffff": "808080") << "\" stroke=\"#ffffff\" />\n";
              }
            }
          }
        }
      }
      p += 2;
    }
  }
}

static void XMLCALL endElement(void* userData, const char* name)
{
  if (!strcmp(name, "g")) {
  	a = d = 1;
  	b = c = e = f = 0;
  	country_cd_1 = country_cd_2 = 0;  	
  }
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
  	perror("Usage:\n\n\tdodecaglobe face mapfile.svg\n\nface is a number between 1 (North pole) and 12 (South pole).");
  	return -1;
  }
  switch (atoi(argv[1])) {
  	case  1: lambda0 =      M_PI / 10; polar = true; face_rotation_angle =      M_PI / 10; break;
  	case  2:                                         face_rotation_angle =  2 * M_PI / 5; break;
  	case  3: lambda0 =  2 * M_PI / 5;                                                     break;
  	case  4: lambda0 =  4 * M_PI / 5;                face_rotation_angle =  4 * M_PI / 5; break;
  	case  5: lambda0 = -4 * M_PI / 5;                face_rotation_angle = -2 * M_PI / 5; break;
  	case  6: lambda0 = -2 * M_PI / 5;                face_rotation_angle = -4 * M_PI / 5; break;
  	case  7: lambda0 =  1 * M_PI / 5; south = true; break;
  	case  8: lambda0 =  3 * M_PI / 5; south = true;  face_rotation_angle =  2 * M_PI / 5; break;
  	case  9: lambda0 =  5 * M_PI / 5; south = true;  face_rotation_angle =  2 * M_PI / 5; break;
  	case 10: lambda0 = -3 * M_PI / 5; south = true;  face_rotation_angle = -2 * M_PI / 5; break;
  	case 11: lambda0 = -1 * M_PI / 5; south = true;  face_rotation_angle =  4 * M_PI / 5; break;
    case 12: polar = true; south = true;             face_rotation_angle = -2 * M_PI / 5; break;
    default:
    	std::cerr << "Expected face between 1 and 12\n";
    	return -2;
  };
  lambda0 += M_PI / 5;
  face_rotation_angle += M_PI / 10;
  std::ifstream f(argv[2], std::ios::in);
  if (!f) {
    std::cerr << "Can't open input file " << argv[1] << "\n";
    return -1;
  }
  std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
  	"<svg version=\"1.1\"\n"
    "  baseprofile=\"full\"\n"
    "  xmlns=\"http://www.w3.org/2000/svg\"\n"
    "  xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
    "  xmlns:ev=\"http://www.w3.org/2001/xml-events\"\n"
    "  width=\"" << printable_width_in_px << "\"\n"
    "  height=\"" << printable_height_in_px << "\"\n"
    "  viewBox=\"0 0 " << printable_width_in_px << " " << printable_height_in_px << "\"\n"
    ">\n";

  // Clip path to avoid printing outside the face
  std::cout << "<defs>\n<clipPath id=\"ClipPath5\">\n<polygon points=\"";
  double angle = 18 * M_PI / 180;
  double angle_1 = angle; 
  for (int i = 0; i < 5; i++) {
  	if (i) {
  	  std::cout << " ";	
  	}
    std::cout << (pentagon_center.real() + pentagon_radius_in_pixels * sin(angle)) << "," 
              << (pentagon_center.imag() + pentagon_radius_in_pixels * cos(angle));
  	angle += 2 * M_PI / 5;
  }
  std::cout << "\" />\n</clipPath>\n</defs>\n";

  // White background
  std::cout << "<rect width=\"" << printable_width_in_px << "\" height=\"" << printable_height_in_px << "\" fill=\"#ffffff\" />\n";

  // Folding lines
  double slope = printable_width_in_px * tan(18 * M_PI / 180) / 2;
  double vertical_shift = (w + s) * output_scaling_factor / 2;
  for (int i = -1; i <= 1; i++) {
    std::cout << "<line x1=\"" << 0 
              << "\" y1=\""    << (printable_height_in_px / 2 + slope - i * vertical_shift) 
              << "\" x2=\""    << printable_width_in_px 
              << "\" y2=\""    << (printable_height_in_px / 2 - slope - i * vertical_shift)
              << "\" stroke=\"black\" stroke-width=\"1\" stroke-dasharray=\"10,10\" />\n";
    slope = -slope;
  }
  const double angle_3 = 18 * M_PI / 180 + 4 * M_PI / 5;
  const double delta_y = (printable_width_in_px - 
    (pentagon_center.real() + pentagon_radius_in_pixels * sin(angle_3))) / tan(36 * M_PI / 180);
  std::cout << "<line x1=\"" << (pentagon_center.real() + pentagon_radius_in_pixels * sin(angle_3)) 
            << "\" y1=\""    << (pentagon_center.imag() + pentagon_radius_in_pixels * cos(angle_3)) 
            << "\" x2=\""    << printable_width_in_px 
            << "\" y2=\""    << ( pentagon_center.imag() + pentagon_radius_in_pixels * cos(angle_3) + delta_y)
            << "\" stroke=\"black\" stroke-width=\"1\" stroke-dasharray=\"10,10\" />\n";
  std::cout << "<line x1=\"" << (mirror_center.real() - pentagon_radius_in_pixels * sin(angle_3)) 
            << "\" y1=\""    << (mirror_center.imag() - pentagon_radius_in_pixels * cos(angle_3)) 
            << "\" x2=\""    << 0 
            << "\" y2=\""    << (mirror_center.imag() - pentagon_radius_in_pixels * cos(angle_3) - delta_y)
            << "\" stroke=\"black\" stroke-width=\"1\" stroke-dasharray=\"10,10\" />\n";
  for (int i = 0; i < 5; i++) {
  	double angle_2 = angle_1 + 2 * M_PI / 5;
  	if (i == 3 || i == 0) {
      std::cout << "<line x1=\"" << (pentagon_center.real() + pentagon_radius_in_pixels * sin(angle_1)) 
                << "\" y1=\""    << (pentagon_center.imag() + pentagon_radius_in_pixels * cos(angle_1)) 
                << "\" x2=\""    << (pentagon_center.real() + pentagon_radius_in_pixels * sin(angle_2)) 
                << "\" y2=\""    << (pentagon_center.imag() + pentagon_radius_in_pixels * cos(angle_2))
                << "\" stroke=\"black\" stroke-width=\"1\" stroke-dasharray=\"10,10\" />\n";
      std::cout << "<line x1=\"" << (mirror_center.real() - pentagon_radius_in_pixels * sin(angle_1)) 
                << "\" y1=\""    << (mirror_center.imag() - pentagon_radius_in_pixels * cos(angle_1)) 
                << "\" x2=\""    << (mirror_center.real() - pentagon_radius_in_pixels * sin(angle_2)) 
                << "\" y2=\""    << (mirror_center.imag() - pentagon_radius_in_pixels * cos(angle_2))
                << "\" stroke=\"black\" stroke-width=\"1\" stroke-dasharray=\"10,10\" />\n";
    }
    angle_1 = angle_2;
  }

  std::cout << "<g clip-path=\"url(#ClipPath5)\">\n"; // Start clip path
  XML_Parser parser = XML_ParserCreate(NULL);
  XML_SetUserData(parser, NULL);
  XML_SetElementHandler(parser, startElement, endElement);
  char buf[BUFSIZ];
  while (!f.eof()) {
  	f.read(buf, sizeof(buf));
  	if (XML_Parse(parser, buf, f.gcount(), f.eof() ? 1 : 0) == XML_STATUS_ERROR) {
  	  std::cerr << XML_ErrorString(XML_GetErrorCode(parser)) 
  	    << " at line " << XML_GetCurrentLineNumber(parser)
  	    << " column " << XML_GetCurrentColumnNumber(parser)
  	    << "\n";
      break;
  	}
  }
  XML_ParserFree(parser);
  std::cout << "</g>\n"; // End clip path

  if (polar && south) { // Patch the hole at the South pole
    std::cout << "<circle r=\"30\" fill=\"#808080\""
                 " cx=\"" << pentagon_center.real() << 
                 "\" cy=\"" << pentagon_center.imag() << "\" />\n";
  }

  std::cout << "</svg>\n";
  return 0;
}
