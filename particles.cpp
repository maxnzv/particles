#include <iostream>
#include <math.h>
#include <time.h>

#include <mutex>

#include <thread>
#include <vector>

#include <xcb/xcb.h>

// Particles amount
const int N = 20;
// Threads amount
const int threads = 2;
std::vector<std::thread> VThread;
// Gravity constant
const double G = 0.0001;
const double minX = 0.0;
const double minY = 0.0;
const double minZ = 0.0;
const double maxX = 1024.0;
const double maxY = 700.0;
const double maxZ = 700.0;
const double MASS = 1.0;
// A force decrement to simulate collisions
const double cf = 0.9;
// Minimal distance to reverce the force between particles
const double minDD = 4.0;
const double minD2 = minDD * minDD * minDD;
// Border width that reverses velocity
const int border = 10;
// Standart delay between display iterations 50 FPS
timespec delay = { 0, 20000000 };
//timespec delay = { 1, 0 };
int finish = 0;
// Calculations per second
long cps = 0;
// Are we ready to move
int movePending;
std::mutex moveMutex;
std::mutex cpsMutex;

class TParticle {
    // Mass
    double m;
    // Coordinates
    double x, y, z;
    // Velocity
    double vx, vy, vz;
    // Force
    double fx, fy, fz;
    // Coordinates for drawing
    xcb_point_t *xcbpoint;
  public:
    TParticle () {
      m = x = y = z = vx = vy = vz = fx = fy = fz = 0;
      xcbpoint = NULL;
    }
    void Init (double X, double Y, double M, xcb_point_t *point) {
      x = X;
      y = Y;
      z = X;
      m = M;
      xcbpoint = point;
      xcbpoint->x = (int)x;
      xcbpoint->y = (int)y;
      vx = vy = vz = fx = fy = fz = 0.0;
    }

    // Save coordinates to future drawing
    void Draw (void) {
      xcbpoint->x = (int)x;
      xcbpoint->y = (int)y;
    }

    // Move the particle according on force, time and initial speed
    // Use clock_gettime with CLOCK_MONOTONIC
    void Move ( long t ) {
      // Some optimization
      double tm = t / m;
      double mt22 = tm*t / 2.0;

      // Border protection
      if ( x < minX+border or x > maxX-border ) vx = -vx;
      if ( y < minY+border or y > maxY-border ) vy = -vy;
      if ( z < minZ+border or z > maxZ-border ) vz = -vz;

      x += vx * t + fx * mt22;
      y += vy * t + fy * mt22;
      z += vz * t + fz * mt22;
      vx += fx * tm;
      vy += fy * tm;
      vz += fz * tm;
      fx = fy = fz = 0.0;
    }

    // Calculate force between two particles and add it to both
    void CalcForce (TParticle *neighbor) {
      double dx = x - neighbor->x;
      double dy = y - neighbor->y;
      double dz = z - neighbor->z;
      double m1 = neighbor->m;
      double d2 = dx*dx + dy*dy + dz*dz; // distance ^ 3 (Pifagor)
      double dd = cbrt(d2);
      double f = G * m*m1/d2;
      if (d2 < minD2) { f = -f*cf; }
      double dfx = f * dx / dd;
      double dfy = f * dy / dd;
      double dfz = f * dz / dd;
      fx -= dfx;
      fy -= dfy;
      fz -= dfz;
      neighbor->fx += dfx;
      neighbor->fy += dfy;
      neighbor->fz += dfz;
    }
};

class TPArray {
    TParticle *container;
    xcb_point_t *xpoints;
  public:
    TPArray (xcb_point_t *xp) {
      xpoints = xp;
      container = new TParticle[N];
      for (int i=0; i<N-1; i++) {
        container[i].Init ((maxX-border*2)*i/N + sqrt(i) + border,
                           (maxY-border*2)*i/N + border, MASS, xpoints + i);
      }
      container[N-1].Init ((maxX-border*2) + border,
                         border * 2, MASS * 2, xpoints + N-1);
    }

    // Calculate forces among "amount" particles original:
    // from start and neighbors from nstart
    void Calculate (int start = 0, int amount = N) {
      for (int i=start; i<start+amount-1; i++) {
        for (int j=i+1; j<N; j++) {
          container[i].CalcForce (container + j);
        }
      }
    }

    // Move every particle according applied forces
    void Move (long t, int start = 0, int amount = N) {
      for (int i=start; i<start+amount; i++) {
        container[i].Move (t);
      }
    }

    // Save every particle position
    void Draw (void) {
      for (int i=0; i<N; i++) {
        container[i].Draw ();
      }
    }
    ~TPArray (){
      delete container;
    }
};

void CalcAndMove (TPArray *ppa, int start = 0, int amount = N) {
  while (!finish) {
    cpsMutex.lock();
    cps++;
    cpsMutex.unlock();
    // Particles in this thread are not ready to move
    moveMutex.lock();
    movePending++;
    moveMutex.unlock();
    ppa->Calculate(start, amount);
    // Particles are ready to move
    moveMutex.lock();
    movePending--;
    moveMutex.unlock();
    // Wait for other threads to complete calculations
    while (movePending) { 
      std::this_thread::yield();
    }
    ppa->Move(1, start, amount);
  }
}

int main () {
  xcb_connection_t    *connection;
  xcb_screen_t        *screen;
  xcb_drawable_t       win;
  xcb_gcontext_t       foreground;
  xcb_gcontext_t       background;
  xcb_generic_event_t *event;
  uint32_t             mask = 0;
  uint32_t             values[2];

  /* geometric objects */
  xcb_point_t          points[N];
  TPArray tpa(points);

  for (int i=0; i<N; i++) {
    std::cout<<points[i].x<<"-"<<points[i].y<<std::endl;
  }

  /* Open the connection to the X server */
  connection = xcb_connect (NULL, NULL);

  /* Get the first screen */
  screen = xcb_setup_roots_iterator (xcb_get_setup (connection)).data;

  /* Create black (foreground) graphic context */
  win = screen->root;

  foreground = xcb_generate_id (connection);
  mask = XCB_GC_FOREGROUND | XCB_GC_GRAPHICS_EXPOSURES;
  values[0] = screen->white_pixel;
  values[1] = 0;
  xcb_create_gc (connection, foreground, win, mask, values);

  background = xcb_generate_id (connection);
  mask = XCB_GC_FOREGROUND | XCB_GC_GRAPHICS_EXPOSURES;
  values[0] = screen->black_pixel;
  values[1] = 0;
  xcb_create_gc (connection, background, win, mask, values);

  /* Ask for our window's Id */
  win = xcb_generate_id(connection);

  /* Create the window */
  mask = XCB_CW_BACK_PIXEL | XCB_CW_EVENT_MASK;
  values[0] = screen->black_pixel;
  values[1] = XCB_EVENT_MASK_EXPOSURE;
  xcb_create_window (connection,                             /* Connection          */
                     XCB_COPY_FROM_PARENT,          /* depth               */
                     win,                           /* window Id           */
                     screen->root,                  /* parent window       */
                     0, 0,                          /* x, y                */
                     int(maxX), int(maxY),                      /* width, height       */
                     10,                            /* border_width        */
                     XCB_WINDOW_CLASS_INPUT_OUTPUT, /* class               */
                     screen->root_visual,           /* visual              */
                     mask, values);                 /* masks */

  // Some magic to close window
  xcb_intern_atom_cookie_t cookie = xcb_intern_atom(connection, 1, 12, "WM_PROTOCOLS");
  xcb_intern_atom_reply_t* reply = xcb_intern_atom_reply(connection, cookie, 0);
  xcb_intern_atom_cookie_t cookie2 = xcb_intern_atom(connection, 0, 16, "WM_DELETE_WINDOW");
  xcb_intern_atom_reply_t* reply2 = xcb_intern_atom_reply(connection, cookie2, 0);
  xcb_change_property(connection, XCB_PROP_MODE_REPLACE, win, (*reply).atom, 4, 32, 1, &(*reply2).atom);

  /* Map the window on the screen */
  xcb_map_window (connection, win);

  /* We flush the request */
  xcb_flush (connection);

  // Creating a mutex and threads
  moveMutex.lock();
  movePending = 0;
  moveMutex.unlock();
  for (int i=0;i<threads;i++) {
    std::cout<<i<<std::endl;
    VThread.push_back(std::thread (CalcAndMove, &tpa, N*i/threads, N/threads));
  }

  //std::thread cm (CalcAndMove, &tpa, 0, N/3);
  //std::thread cm1 (CalcAndMove, &tpa, N/3+1, N/3);
  //std::thread cm2 (CalcAndMove, &tpa, N/3*2+2, N/3);

  int n = 0;

  while (1) {
    if ((event = xcb_poll_for_event (connection))== NULL) {
      clock_nanosleep (CLOCK_MONOTONIC, 0, &delay, &delay);

      //tpa.Draw();
      //for (int i=0; i<N; i++) {
        //std::cout<<points[i].x<<std::endl;
      //}

      xcb_poly_point (connection, XCB_COORD_MODE_ORIGIN, win, background, N, points);
      tpa.Draw();
      xcb_poly_point (connection, XCB_COORD_MODE_ORIGIN, win, foreground, N, points);
      xcb_flush (connection);
      if (n == 50) {
        moveMutex.lock();
        std::cout<<"tick: "<<cps<<", "<<movePending<<std::endl;
        moveMutex.unlock();
        n = 0;
	cpsMutex.lock();
        cps = 0;
	cpsMutex.unlock();
        //for (int i=0; i<N; i++) {
          //std::cout<<points[i].x<<"-"<<points[i].y<<std::endl;
        //}
      }
      else {n++;}
      continue;
    }
    switch (event->response_type & ~0x80) {
      case XCB_EXPOSE: {
        /* We draw the points */
        xcb_poly_point (connection, XCB_COORD_MODE_ORIGIN, win, foreground, N, points);

        /* We flush the request */
        xcb_flush (connection);

        break;
      }
      case XCB_CLIENT_MESSAGE: {
          // The magic continues
          if((*(xcb_client_message_event_t*)event).data.data32[0] == (*reply2).atom) {
            finish = 1;
            for (auto& th: VThread) {
              th.join();
            }
            //cm.join();
            //cm1.join();
            //cm2.join();
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
    free (event);
  }

  return 0;
}
