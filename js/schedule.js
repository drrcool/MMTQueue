$(document).ready(function() {
// Page is ready, initialize the calendar
$('#calendar').fullCalendar({
    events: 'schedule.json',
    header: {
        left: ',prev,next today',
        center: 'title',
        right: 'month,agendaWeek,agendaDay'
    },
    defaultView: 'month',
    allDaySlot: false,
    minTime: "00:00:00",
    maxTime: "16:00:00",
    eventTextColor: "black",
    displayEventTime: true,

    eventRender: function(event, element) {
        element.find('.fc-title').html('<center>'+
            element.find('.fc-title').text() +
            '</center>');
    },

    dayClick: function(date, jsEvent, view) {
        $('#calendar').fullCalendar('select', date);
        $('#calendar').fullCalendar('gotoDate', date);
        $('#calendar').fullCalendar('changeView', 'agendaWeek');
    },

    eventClick: function(event, jsEvent, view) {
        $('#modalTitle').html(event.title);
        $('#modal_ra').html(event.ra);
        $('#modal_obstype').html(event.obstype);
        $('#modal_pi').html(event.PI);
        $('#modal_dec').html(event.dec);
        $('#modal_mag').html(event.mag);
        $('#modal_mask').html(event.mask);
        $('#modal_filter').html(event.filter);
        $("#modal_grism").html(event.grism);
        $('#modal_gain').html(event.gain);
        $('#modal_readtab').html(event.readtab);
        $('#modal_exptime').html(event.exptime*60.);
        $('#modal_nexp').html(event.repeats);
        $('#modal_nvisits_scheduled').html(event.n_visits_scheduled);
        $('#modal_dithersize').html(event.dithersize);
        $("#modal_nexposures").html(event.n_visits_scheduled*event.repeats);
        $('#modal_totalexptime').html(event.n_visits_scheduled*event.repeats*event.exptime*60);
        $('#modal_photometric').html(event.photometric);
        $('#modal_seeing').html(event.seeing);
        console.log(event);
        $('#fullCalModal').modal();
        return false;
    }
});
});
