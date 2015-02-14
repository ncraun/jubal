/*
    Copyright (C) 2014-2015 Nate Craun

    This file is part of Jubal.

    Jubal is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jubal is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jubal.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef JUBAL_UNDO_H
#define JUBAL_UNDO_H

/*
	Undo functionality. Each possible action in the program is represented
	by an action. These actions are placed on the stack whose size is
	limited by a possibly user defined amount. When the users selects undo
	the most recent action is popped off the stack, and a reverse action 
	is performed. The action is then placed on a redo stack, which will be
	popped when the user "redoes" actions. The redo stack is the same size
	as the undo stack. When the stacks reach their maximum sizes, the 
	oldest item on the stack will be removed, and the new action takes its
	place.
*/

/*
	Actions:

	Note Edit
	Transpose
	Delete
	Insert
	Shift
	Cut
	Param Change

	Node Delete
	Node Add
	Node Change Def
	-- Note: To support these changes, we can't free nodes when the user
		deletes them. We must wait for the undo to go out of scope.
*/

#endif
