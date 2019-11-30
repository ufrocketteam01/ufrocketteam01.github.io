jQuery(document).ready(function(){
    "use strict"
    $("a.smooth-scroll").click(function (event) {

		event.preventDefault();

		var section = $(this).attr("href");

		$('html, body').animate({
		  scrollTop: $(section).offset().top - -2
		}, 1250, "easeInOutExpo");
  	});
   	new WOW().init();
});
$("#team-members").owlCarousel({
	items:4,
	autoplay:true,
	smartSpeed:100,
	loop:true,
	autoplayHoverPause:true,
	 
});
// $('.counter').counterUp({
// 		  delay: 10,
// 		  time: 4000
// 	  });


$(document).ready(function () {

    $('[data-toggle="popover"]').popover() });
 


$('#part').popover('update')
var clickTop,clickLeft=0;
    $(document).click(function(e){
        clickTop =e.pageY;
        clickLeft =e.pageX;

    });
    $().ready(function(){
        var popovers=$('[data-toggle="popover"]');

        popovers.popover({
            placement: 'bottom center',
            html:true,
            trigger:'focus'
        }).on("shown.bs.popover", function(e){
            $('.popover').css({top:clickTop-100,left:clickLeft-130});
        })

    });

$('html').on('click', function(e) {
  if (typeof $(e.target).data('original-title') == 'undefined' &&
     !$(e.target).parents().is('.popover.in')) {
    $('[data-original-title]').popover('hide');
  }
});
$('html').popover({
                selector: '[rel=popover]',
                trigger: "click"
            }).on("show.bs.popover", function(e){
                // hide all other popovers
                $("[rel=popover]").not(e.target).popover("destroy");
                $(".popover").remove();                    
            });







