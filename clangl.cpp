#include <ClanLib/application.h>
#include <ClanLib/core.h>
#include <ClanLib/sound.h>
#include <ClanLib/gui.h>
#include <ClanLib/display.h>
#include <ClanLib/swrender.h>
#include <ClanLib/gl.h>

#include <math.h>


class PushButton : public CL_Window
{
public:
	PushButton(CL_GUIManager &manager, CL_ResourceManager &application_resources) :
	        	CL_Window(&manager, CL_GUITopLevelDescription("PushButton", CL_Rect(256 + 16, 256 + 16, CL_Size(256, 180)), false))
	        {
	        	CL_GraphicContext gc = get_gc();
	        	test_image = CL_Image(gc, "tux", &application_resources);

	        	set_draggable(true);

	        	CL_Rect client_area = get_client_area();

	        	pushbutton1 = new CL_PushButton(this);
	        	pushbutton1->set_geometry(CL_Rect(client_area.left + 11, client_area.top + 10, CL_Size(128, 40)));
	        	pushbutton1->set_text("Push Button");
	        	pushbutton1->set_flat(false);
	        	pushbutton1->func_clicked().set(this, &PushButton::on_clicked, pushbutton1);

	        	int label_xpos = client_area.left + 31;
	        	int yoffset = client_area.top + 80;
	        	CL_Size label_size(50, 15);
	        	const int gap = 16;

	        	checkbox_flat = new CL_CheckBox(this);
	        	checkbox_flat->set_geometry(CL_Rect(client_area.left + 11, yoffset, CL_Size(100, 15)));
	        	checkbox_flat->func_checked().set(this, &PushButton::on_checked_flat, checkbox_flat);
	        	checkbox_flat->func_unchecked().set(this, &PushButton::on_unchecked_flat, checkbox_flat);
	        	checkbox_flat->set_text("Flat");

	        	yoffset+=gap;
	        	checkbox_icon = new CL_CheckBox(this);
	        	checkbox_icon->set_geometry(CL_Rect(client_area.left + 11, yoffset, CL_Size(100, 15)));
	        	checkbox_icon->func_checked().set(this, &PushButton::on_checked_icon, checkbox_icon);
	        	checkbox_icon->func_unchecked().set(this, &PushButton::on_unchecked_icon, checkbox_icon);
	        	checkbox_icon->set_text("Icon");

	        	yoffset+=gap;
	        	checkbox_toggle = new CL_CheckBox(this);
	        	checkbox_toggle->set_geometry(CL_Rect(client_area.left + 11, yoffset, CL_Size(100, 15)));
	        	checkbox_toggle->func_checked().set(this, &PushButton::on_checked_toggle, checkbox_toggle);
	        	checkbox_toggle->func_unchecked().set(this, &PushButton::on_unchecked_toggle, checkbox_toggle);
	        	checkbox_toggle->set_text("Enable Toggle");

	        	yoffset+=gap;
	        	checkbox_disable = new CL_CheckBox(this);
	        	checkbox_disable->set_geometry(CL_Rect(client_area.left + 11, yoffset, CL_Size(100, 15)));
	        	checkbox_disable->func_checked().set(this, &PushButton::on_checked_disable, checkbox_disable);
	        	checkbox_disable->func_unchecked().set(this, &PushButton::on_unchecked_disable, checkbox_disable);
	        	checkbox_disable->set_text("Disable");

	        	int xoffset = client_area.left + 36;
	        	yoffset = client_area.top + 60;

//	        	info_clicked = new Info(this);
//	        	info_clicked->set(xoffset, yoffset, "Clicked");
	        }

        virtual ~PushButton() {}

private:




        void on_clicked(CL_PushButton *pushbutton)
        {
//        	info_clicked->activate();
        }

        void on_checked_disable(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_enabled(false);
        }

        void on_unchecked_disable(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_enabled(true);
        }

        void on_checked_icon(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_icon(test_image);
        }

        void on_unchecked_icon(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_icon(CL_Image());
        }

        void on_checked_toggle(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_toggle(true);
        }

        void on_unchecked_toggle(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_toggle(false);
        }

        void on_checked_flat(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_flat(false);
        }

        void on_unchecked_flat(CL_CheckBox *checkbox)
        {
        	pushbutton1->set_flat(true);
        }


private:
        CL_PushButton *pushbutton1;
        CL_CheckBox *checkbox_disable;
        CL_CheckBox *checkbox_icon;
        CL_CheckBox *checkbox_flat;
        CL_CheckBox *checkbox_toggle;

        CL_Image test_image;

//        Info *info_clicked;
};






class GUI
{
public:
	GUI(CL_DisplayWindow * wnd)
	{


		//resources_internal = CL_ResourceManager("../CommonCode/Resources/resources.xml");

		CL_GraphicContext gc = wnd->get_gc();
		fps_font = CL_Font(gc, "Tahoma", 24);







		//std::cout << gui_manager.impl.get() << "\n";
		gui_manager.set_window_manager(wm);

		// Use a texture group to store all the gui textures

		CL_TextureGroup texture_group(CL_Size(1024, 1024));
		wm.set_texture_group(texture_group);    // Note: This line is optional

		resources_gui = CL_ResourceManager(gui->get_resources_location());

		theme.set_resources(resources_gui);
		gui_manager.set_theme(theme);

		gui_manager.set_css_document(gui->get_theme_location());

		// Since this example rebuilds the gui when the theme changes, we have to manually delete the components
		// (We could have instead recreated the CL_GUIManager, so it's destructor deletes them)

		pushbutton.reset(new PushButton(gui->get_gui_manager(), gui->get_resources_internal()));
	}


	const char *get_theme_location()
	{

		return "Resources/GUIThemeLuna/theme.css";
	}

	const char *get_resources_location()
	{

		return "Resources/GUIThemeLuna/resources.xml";
	}


	bool run(CL_GraphicContext &gc)
	{
			wm.process();
			wm.draw_windows(gc);

			return true;
	}

    bool run( CL_DisplayWindow *wnd )
    {
    	static int total_time = 0, fps_count = 0, last_fps= 0;
    	static int start_time = 0;

    	if (start_time == 0)
    	{
    		start_time = CL_System::get_time();
    	}

    	CL_GraphicContext gc = wnd->get_gc();

    	gc.set_map_mode(CL_MapMode(cl_map_2d_upper_left));

    	CL_Draw::gradient_fill(gc, wnd->get_viewport(), CL_Gradient(CL_Colorf(0.4f, 0.4f, 0.4f, 1.0f), CL_Colorf(0.0f, 0.0f, 0.0f, 1.0f)));

    	run_manager(gc);

    	CL_String fps = cl_format("FPS: %1", last_fps);
    	fps_font.draw_text(gc, gc.get_width() - 100 - 2, 24 - 2, fps, CL_Colorf(0.0f, 0.0f, 0.0f, 1.0f));
    	fps_font.draw_text(gc, gc.get_width() - 100, 24, fps, CL_Colorf::white);

    	fps_font.draw_text(gc, 24, gc.get_height() - 16, "Rendering GUI onto a CL_Texture, then onto the display window.", CL_Colorf(1.0f, 1.0f, 1.0f, 1.0f));

    	fps_count++;
    	int time = CL_System::get_time();
    	total_time += time - start_time;
    	start_time = time;
    	if(total_time >= 1000)
    	{
    		last_fps = fps_count;
    		total_time -= 1000;
    		fps_count = 0;
    	}




    	return true;
    }


        CL_DisplayWindow *get_window();
        //App *get_app() { return app; }
        CL_GUIManager &get_gui_manager() { return gui_manager; }
        CL_ResourceManager &get_resources_internal() { return resources_internal; }

        //Theme::gui_theme get_theme() {return current_theme;}



private:




        void run_manager(CL_GraphicContext &gc)
        {



			run(gc);


        }



//        void gui_repaint();
//        void gui_exec();
//        void run_manager(CL_GraphicContext &gc);
//        void reset_manager();

private:
        CL_GUIManager gui_manager;
        CL_ResourceManager resources_internal;
        //App *app;
        CL_Font fps_font;


        CL_ResourceManager resources_gui;
		CL_GUIThemeDefault theme;
		GUI *gui;
		CL_DisplayWindow *window_ptr;
		CL_GUIWindowManagerTexture wm;


		CL_UniquePtr<PushButton> pushbutton;

        //Theme::gui_theme new_theme;
        //Theme::gui_theme current_theme;

        //Balls balls;

};




struct cube
{
	struct
	{
		float pos[3];
		float col[3];
	}ver[8];

	struct
	{
		unsigned int ver[4];
	} quad[6];


	cube() {


		float vertexPosDat[8][3]=
		{
				{-10, 10, 10}, //left,top,front
				{10,  10, 10}, //right,top,front
				{10,  10,-10}, //right,top,back
				{-10, 10,-10}, //left, top,back
				{-10,-10, 10}, //left,bottom,front
				{10, -10, 10}, //right,bottom,front
				{10, -10,-10}, //right,bottom,back
				{-10,-10,-10}  //left,bottom,back
		};

		//defines the colour of each vertex
		float vertexColDat[8][3]=
		{
				{0.5,0  ,0  }, //dark red
				{1  ,1  ,0.3}, //yellow
				{1  ,0  ,0  }, //red
				{0.5,1  ,0.2}, //dull yellow??
				{1  ,1  ,0  }, //yellow
				{0.9,0.5,0  }, //orange
				{1  ,0.9,0.1}, //yellow
				{1  ,0  ,0  }, //red
		};

		//defines the vertexes of each quad in anti-clockwise order
		unsigned int quadVerDat[6][4]=
		{
				{0,1,2,3}, //top
				{0,3,7,4}, //left
				{3,2,6,7}, //back
				{2,1,5,6}, //right
				{0,4,5,1}, //front
				{4,7,6,5}, //bottom
		};


		int a,b;
		//put the vertex data into the cube struct
		for (a=0;a<8;++a)
		{
			for (b=0;b<3;++b)
			{
				ver[a].pos[b]=vertexPosDat[a][b];
				ver[a].col[b]=vertexColDat[a][b];
			}
		}

		//put the quad data into the cube struct
		for (a=0;a<6;++a)
		{
			for (b=0;b<4;++b)
			{
				quad[a].ver[b]=quadVerDat[a][b];
			}
		}
	}

	void redraw() {


		unsigned int currentVer;



		for (int a=0;a<6;++a)
		{
			for ( int b=0;b<4;++b)
			{
				currentVer=quad[a].ver[b];

				glColor3fv(ver[ currentVer ].col);
				glVertex3fv(ver[ currentVer ].pos);
			}
		}

	}

};



class test {


public:
	static void glhPerspectivef2(float *matrix, float fovyInDegrees, float aspectRatio,
			float znear, float zfar)
	{
		float ymax, xmax;
		float temp, temp2, temp3, temp4;
		ymax = znear * tanf(fovyInDegrees * M_PI / 360.0);
		//ymin = -ymax;
		//xmin = -ymax * aspectRatio;
		xmax = ymax * aspectRatio;
		glhFrustumf2(matrix, -xmax, xmax, -ymax, ymax, znear, zfar);
	}
	static void glhFrustumf2(float *matrix, float left, float right, float bottom, float top,
			float znear, float zfar)
	{
		float temp, temp2, temp3, temp4;
		temp = 2.0 * znear;
		temp2 = right - left;
		temp3 = top - bottom;
		temp4 = zfar - znear;
		matrix[0] = temp / temp2;
		matrix[1] = 0.0;
		matrix[2] = 0.0;
		matrix[3] = 0.0;
		matrix[4] = 0.0;
		matrix[5] = temp / temp3;
		matrix[6] = 0.0;
		matrix[7] = 0.0;
		matrix[8] = (right + left) / temp2;
		matrix[9] = (top + bottom) / temp3;
		matrix[10] = (-zfar - znear) / temp4;
		matrix[11] = -1.0;
		matrix[12] = 0.0;
		matrix[13] = 0.0;
		matrix[14] = (-temp * zfar) / temp4;
		matrix[15] = 0.0;
	}


	test() {

	}

	void mainloop() {
		CL_OpenGLWindowDescription desc;
		desc.set_size( CL_Size( 1024, 768 ), true );

		CL_DisplayWindow wnd(desc);

		CL_GraphicContext gc = wnd.get_gc();





		glMatrixMode(GL_PROJECTION);						//hello
		//		gluPerspective(45, //view angle
		//					1.0,	//aspect ratio
		//					10.0, //near clip
		//					200.0);//far clip

		GLfloat matrix[16];
		glhPerspectivef2(matrix, 45, 1.0, 10.0, 200.0 );

		glLoadMatrixf( matrix );

		glMatrixMode(GL_MODELVIEW);




		cube c;
		float rotateBy;

		CL_Font font(gc, "Helvetica", 24);


		GUI gui(&wnd);

		while( true ) {



			rotateBy += 1;
			glEnable(GL_CULL_FACE);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


			glPushMatrix();


			glTranslatef(0,0,-50);
			glRotatef(rotateBy,1,1,0);

			glBegin(GL_QUADS);

			c.redraw();

			glEnd();
			glPopMatrix();


			glDisable( GL_CULL_FACE );
			CL_Draw::line(gc, 0, 0, 100, 100, CL_Colorf( 1.0f ,1.0f ,1.0f ));
			font.draw_text( gc, 100, 100, "bla bla bla", CL_Colorf( 1.0f, 1.0f, 1.0f ));

			gui.run(&wnd);
			wnd.flip();

			CL_System::sleep( 10 );
			CL_KeepAlive::process();
		}

	}



private:
	CL_SetupCore setup_core_;
	CL_SetupSound setup_sound_;



	CL_SetupDisplay display_;
	//CL_SetupSWRender swrender;

	CL_SetupGL setup_gl_;



};


class app {
public:
	static int main(const std::vector<CL_String> &args) {

		test t1;
		t1.mainloop();

		return 0;
	}
};

CL_ClanApplication app(&app::main);
