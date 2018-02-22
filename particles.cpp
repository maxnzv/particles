#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>

#include <xcb/xcb.h>

// Particles amount
const int N = 16;
// Gravity constant
const double G = 10.0;
const double minX = 10.0;
const double minY = 10.0;
const double maxX = 700.0;
const double maxY = 700.0;
// Minimal distance to reverce the force between particles
const double minD2 = 5.0 * 5.0;
// Border width that reverses velocity
const int border = 10;
// Standart delay between display iterations
timespec delay = { 1, 0 };

class TParticle {
    // Mass
    double m;
    // Coordinates
    double x, y;
    // Velocity
    double vx, vy;
    // Force
    double fx, fy;
    // Coordinates for drawing
    xcb_point_t *xcbpoint;
  public:
    TParticle () {
      m = x = y = vx = vy = fx = fy = 0;
      xcbpoint = NULL;
    }
    void Init (double X, double Y, double M, xcb_point_t *point) {
      x = X;
      y = Y;
      m = M;
      xcbpoint = point;
      xcbpoint->x = (int)x;
      xcbpoint->y = (int)y;
      vx = vy = fx = fy = 0.0;
    }

    double GetX (void) { return x; }
    double GetY (void) { return y; }
    double GetM (void) { return m; }

    // Move the particle according on force, time and initial speed
    // Use clock_gettime with CLOCK_MONOTONIC
    void Move ( long t ) {
      // Some optimization
      double tm = t / m;
      double mt22 = tm*t / 2.0;
      x += vx * t + fx * mt22;
      y += vy * t + fy * mt22;
      vx += fx * tm;
      vy += fy * tm;
      fx = fy = 0.0;
      // Border protection
      if ( x < minX or x > maxX ) vx = -vx;
      if ( y < minY or y > maxY ) vy = -vy;
      xcbpoint->x = (int)x;
      xcbpoint->y = (int)y;
    }

    // Calculate force between two particles and add it to both
    void CalcForce (TParticle *neighbor) {
      double dx = x - neighbor->GetX();
      double dy = y - neighbor->GetY();
      double m1 = neighbor->GetM();
      double d2 = dx*dx + dy*dy; // distance ^ 2 (Pifagor)
      double f = G * m*m1/d2;
      if (d2 < minD2) { f = -f; }
      double dfx = f * dx / dy;
      double dfy = f * dy / dx;
      fx -= dfx;
      fy -= dfy;
      neighbor->fx += dfx;
      neighbor->fy += dfy;
    }
};

class TPArray {
    TParticle *container;
    xcb_point_t *xpoints;
  public:
    TPArray (xcb_point_t *xp) {
      xpoints = xp;
      container = new TParticle[N];
      for (int i=0; i<N; i++) {
        container[i].Init (maxX*i/N, maxY*i/N, (double)i, xpoints + i);
      }
    }
    void Calculate (void) {
      for (int i=0; i<N-1; i++) {
        for (int j=i; j<N; j++) {
          container[i].CalcForce (container + j);
	}
      }
    }
    void Move (long t) {
      for (int i=0; i<N-1; i++) {
        container[i].Move (t);
      }
    }
};

int main () {
  xcb_connection_t    *c;
  xcb_screen_t        *screen;
  xcb_drawable_t       win;
  xcb_gcontext_t       foreground;
  xcb_generic_event_t *e;
  uint32_t             mask = 0;
  uint32_t             values[2];

  /* geometric objects */
  xcb_point_t          points[N];
  TPArray tpa(points);

  for (int i=0; i<N; i++) {
    std::cout<<points[i].x<<std::endl;
  }

  /* Open the connection to the X server */
  c = xcb_connect (NULL, NULL);

  /* Get the first screen */
  screen = xcb_setup_roots_iterator (xcb_get_setup (c)).data;

  /* Create black (foreground) graphic context */
  win = screen->root;

  foreground = xcb_generate_id (c);
  mask = XCB_GC_FOREGROUND | XCB_GC_GRAPHICS_EXPOSURES;
  values[0] = screen->white_pixel;
  values[1] = 0;
  xcb_create_gc (c, foreground, win, mask, values);

  /* Ask for our window's Id */
  win = xcb_generate_id(c);

  /* Create the window */
  mask = XCB_CW_BACK_PIXEL | XCB_CW_EVENT_MASK;
  values[0] = screen->black_pixel;
  values[1] = XCB_EVENT_MASK_EXPOSURE;
  xcb_create_window (c,                             /* Connection          */
                     XCB_COPY_FROM_PARENT,          /* depth               */
                     win,                           /* window Id           */
                     screen->root,                  /* parent window       */
                     0, 0,                          /* x, y                */
                     int(maxX)+border, int(maxY)+border,                      /* width, height       */
                     10,                            /* border_width        */
                     XCB_WINDOW_CLASS_INPUT_OUTPUT, /* class               */
                     screen->root_visual,           /* visual              */
                     mask, values);                 /* masks */

  // Some magic to close window
  xcb_intern_atom_cookie_t cookie = xcb_intern_atom(c, 1, 12, "WM_PROTOCOLS");
  xcb_intern_atom_reply_t* reply = xcb_intern_atom_reply(c, cookie, 0);
  xcb_intern_atom_cookie_t cookie2 = xcb_intern_atom(c, 0, 16, "WM_DELETE_WINDOW");
  xcb_intern_atom_reply_t* reply2 = xcb_intern_atom_reply(c, cookie2, 0);
  xcb_change_property(c, XCB_PROP_MODE_REPLACE, win, (*reply).atom, 4, 32, 1, &(*reply2).atom);

  /* Map the window on the screen */
  xcb_map_window (c, win);


  /* We flush the request */
  xcb_flush (c);

  while (1) {
    if ((e = xcb_poll_for_event (c))== NULL) {
        clock_nanosleep (CLOCK_MONOTONIC, 0, &delay, &delay);
	std::cout<<"tick"<<std::endl;
        continue;
    }
    switch (e->response_type & ~0x80) {
      case XCB_EXPOSE: {
        /* We draw the points */
        xcb_poly_point (c, XCB_COORD_MODE_ORIGIN, win, foreground, N, points);

        /* We flush the request */
        xcb_flush (c);

        break;
      }
      case XCB_CLIENT_MESSAGE: {
          // The magic continues
          if((*(xcb_client_message_event_t*)e).data.data32[0] == (*reply2).atom) {
              return 0;
	  }
	  break;
      }
      default: {
        /* Unknown event type, ignore it */
        break;
      }
    }
    /* Free the Generic Event */
    free (e);
  }

  return 0;
}
