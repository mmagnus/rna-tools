   function del(t){
      if (window.confirm("Do you really want to delete file: " + t)) {
        $.ajax({
            url: "/del/" + t,
	})

      }
     };

   function rm(){
        $.ajax({
            url: "/run/clear/{{ j.job_id }}",
	})
    };

